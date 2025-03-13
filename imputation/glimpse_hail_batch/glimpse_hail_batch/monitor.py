import asyncio
import argparse
import statistics
from functools import partial
import re

from rich.live import Live
from rich.table import Table
from rich.console import Group
from rich.panel import Panel

import hailtop.batch_client.aioclient as bc
from hailtop.utils import bounded_gather
from hailtop.utils.time import parse_timestamp_msecs, time_msecs, humanize_timedelta_msecs, time_msecs_str


def get_sample_group_id(name: str) -> int:
    return int(name.split('/')[0].split('-')[2])


async def get_true_job_duration_mins(b: bc.Batch, job: dict):
    name = job['name']
    job = await b.get_job(job['job_id'])
    attempts = await job.attempts()
    duration_ms = 0
    start_time = None
    end_time = None
    if not attempts:
        return (name, None, None, None)

    for attempt in attempts:
        duration_ms += attempt['duration_ms']
        _start_time = parse_timestamp_msecs(attempt.get('start_time'))
        start_time = min(start_time, _start_time) if start_time else _start_time
        _end_time = parse_timestamp_msecs(attempt.get('end_time'))
        end_time = min(end_time, _end_time) if end_time else _end_time
    return (name, start_time, end_time, duration_ms / 1000 / 60)


def markup_cell(value: str, state: str, started: bool) -> str:
    if state in ('Success', 'success'):
        return f'[green]{value}[/]'
    if state in ('failure', 'cancelled', 'Failed', 'Cancelled', 'Errored'):
        return f'[red]{value}[/]'
    if started and state in ('running', 'Ready', 'Creating', 'Running'):
        return f'[blue]{value}[/]'
    else:
        return f'[black]{value}[/]'


async def generate_sample_group_table(b: bc.Batch) -> Table:
    table = Table()

    table.add_column("Sample Group ID")
    table.add_column("Job Group ID")
    table.add_column("Start Time")
    table.add_column("End Time")
    table.add_column("Duration (mins)")
    table.add_column("State")
    table.add_column("Number of Jobs")
    table.add_column("Number of Completed Jobs")
    table.add_column("Percent Completed")
    table.add_column("Number of Succeeded Jobs")
    table.add_column("Number of Failed Jobs")
    table.add_column("Number of Cancelled Jobs")
    table.add_column("Average Phasing Job Runtime (mins)")
    table.add_column("Minimum Phasing Job Runtime (mins)")
    table.add_column("Maximum Phasing Job Runtime (mins)")
    table.add_column("Phase Cost")
    table.add_column("Ligate Cost")
    table.add_column("Other Costs")
    table.add_column("Total Cost")
    table.add_column("Total Cost per Sample")

    job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]
    job_groups = [(name, jg) for name, jg in job_groups if name.startswith('sample-group')]
    job_groups.sort(key=lambda x: get_sample_group_id(x[0]))

    for jg_name, jg in job_groups:
        sample_group_id = get_sample_group_id(jg_name)

        status = await jg.status()

        n_samples = int((await jg.attributes())['N'])

        n_jobs = status['n_jobs']
        n_succeeded = status['n_succeeded']
        n_failed = status['n_failed']
        n_cancelled = status['n_cancelled']
        n_completed = status['n_completed']
        percent_completed = 100 * (n_completed / n_jobs)

        total_cost = status['cost']

        child_job_groups = {(await child_jg.attributes()).get('name', ''): child_jg async for child_jg in jg.job_groups()}
        child_job_groups = {name.split('/')[1]: child_jg for name, child_jg in child_job_groups.items()}

        phase_job_group = child_job_groups['phase']

        phase_status = await phase_job_group.status()
        ligate_status = await child_job_groups['ligate'].status()

        phase_cost = phase_status['cost']
        ligate_cost = ligate_status['cost']

        other_costs = max(0, total_cost - phase_cost - ligate_cost)

        all_jobs = [job async for job in jg.jobs(recursive=True)]
        info = await bounded_gather(*[partial(get_true_job_duration_mins, b, job) for job in all_jobs],
                                    cancel_on_error=True)

        start_times = [x[1] for x in info if x[1] is not None]
        end_times = [x[2] for x in info if x[2] is not None]

        phase_duration_min = [x[3] for x in info if 'phase' in x[0] and x[3] is not None]

        start_time_str = None
        end_time_str = None

        if start_times:
            start_time = min(start_times)
            start_time_str = time_msecs_str(start_time)
            if end_times:
                end_time = max(end_times)
                end_time_str = time_msecs_str(end_time)
                job_group_duration = humanize_timedelta_msecs(end_time - start_time)
            else:
                now = time_msecs()
                job_group_duration = humanize_timedelta_msecs(now - start_time)
        else:
            job_group_duration = None

        if phase_duration_min:
            mean_duration_min = round(statistics.mean(phase_duration_min), 1)
            min_duration_min = round(min(phase_duration_min), 1)
            max_duration_min = round(max(phase_duration_min), 1)
        else:
            mean_duration_min = None
            min_duration_min = None
            max_duration_min = None

        state = status['state']
        if state == 'running' and not start_times:
            state = 'pending'

        row = [f"{sample_group_id}", f"{jg.job_group_id}", f"{start_time_str}", f"{end_time_str}", f"{job_group_duration}",  f"{state}", f"{n_jobs}", f"{n_completed}",
               f"{percent_completed:.2f}", f"{n_succeeded}", f"{n_failed}", f"{n_cancelled}", f"{mean_duration_min}",
               f"{min_duration_min}", f"{max_duration_min}",
               f"${phase_cost:.4f}", f"${ligate_cost:.4f}", f"${other_costs:.4f}", f"${total_cost:.4f}", f"${total_cost / n_samples:.4f}"]

        row = [markup_cell(c, status['state'], (n_succeeded + n_failed + n_cancelled) > 0) for c in row]

        table.add_row(*row)

    return table


async def generate_union_table(b: bc.Batch) -> Table:
    table = Table()

    table.add_column("Contig")
    table.add_column("State")

    job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]

    job_groups = [jg for name, jg in job_groups if name == 'union-sample-groups']
    if not job_groups:
        return table

    contig_job_groups = [((await child_jg.attributes())['contig'], child_jg) async for child_jg in job_groups[0].job_groups()]
    contig_job_groups.sort(key=lambda x: int(re.sub('[^0-9]','', x[0])))

    for contig, contig_jg in contig_job_groups:
        status = await contig_jg.status()

        total_cost = status['cost']

        n_succeeded = status['n_succeeded']
        n_failed = status['n_failed']
        n_cancelled = status['n_cancelled']

        jobs = [(await b.get_job(job['job_id'])) async for job in contig_jg.jobs()]
        state = (await jobs[0].status())['state']

        row = [contig, state, f'${total_cost:.2f}']

        row = [markup_cell(c, state, (n_succeeded + n_failed + n_cancelled) > 0) for c in row]

        table.add_row(*row)

    return table


async def get_total_costs(b: bc.Batch) -> Table:
    table = Table()

    table.add_column("Total Cost")
    table.add_column("GLIMPSE Costs")
    table.add_column("Union + Other Costs")

    job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]

    glimpse_cost = sum([(await jg.status())['cost'] for name, jg in job_groups if name.startswith('sample-group')])
    total_cost = (await b.status())['cost']
    union_cost = total_cost - glimpse_cost

    row = [f'${total_cost:.2f}', f'${glimpse_cost:.2f}', f'${union_cost:.2f}']

    table.add_row(*row)

    return table


async def main(batch_id: int):
    batch_client = await bc.BatchClient.create('dummy')
    try:
        b = await batch_client.get_batch(batch_id)
        with Live(refresh_per_second=4) as live:
            while True:
                sample_group_table = await generate_sample_group_table(b)
                union_sample_groups_table = await generate_union_table(b)
                total_costs = await get_total_costs(b)

                panel_group = Group(
                    Panel(sample_group_table,  title='Sample Groups'),
                    Panel(union_sample_groups_table, title='Union Data By Contig'),
                    Panel(total_costs, title='Total Costs')
                )

                live.update(panel_group)
                await asyncio.sleep(15)
    finally:
        await batch_client.close()



if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--batch-id', type=int, required=True)
    args = p.parse_args()

    asyncio.run(main(args.batch_id))
