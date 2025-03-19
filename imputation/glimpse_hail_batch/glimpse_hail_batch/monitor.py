import asyncio
import argparse
import statistics
from functools import partial
import re
from typing import List

from rich.live import Live
from rich.table import Table
from rich.console import Group
from rich.panel import Panel

import hailtop.batch_client.aioclient as bc
from hailtop.utils import bounded_gather
from hailtop.utils.time import parse_timestamp_msecs, time_msecs, humanize_timedelta_msecs, time_msecs_str


SAMPLE_GROUPS = []


def markup_cell(value: str, state: str, started: bool) -> str:
    if state in ('Success', 'success'):
        return f'[green]{value}[/]'
    if state in ('failure', 'cancelled', 'Failed', 'Cancelled', 'Errored'):
        return f'[red]{value}[/]'
    if started and state in ('running', 'Ready', 'Creating', 'Running'):
        return f'[blue]{value}[/]'
    else:
        return f'[black]{value}[/]'


class SampleGroupProgress:
    @staticmethod
    async def initialize(b: bc.Batch, job_group: bc.JobGroup):
        jobs = [job async for job in job_group.jobs(recursive=True)]
        jobs = [bc.Job.submitted_job(b, j['job_id'], _status=None) for j in jobs]
        progress = SampleGroupProgress(job_group, jobs)
        return progress

    def __init__(self, job_group: bc.JobGroup, jobs: List[bc.Job]):
        self.job_group = job_group
        self.jobs = jobs
        self._job_group_id = job_group.job_group_id
        self._sample_group_id = None
        self._status = None
        self._attempts = None
        self._attributes = None
        self._state = None
        self._start_time = None
        self._end_time = None
        self._cost = None
        self._n_jobs = None
        self._n_succeeded = None
        self._n_failed = None
        self._n_cancelled = None
        self._n_completed = None
        self._refresh_task = None
        self._needs_refresh = True

    async def refresh(self):
        if not self._needs_refresh:
            return

        self._status = await self.job_group.status()
        self._attributes = await self.job_group.attributes()
        self._sample_group_id = int(self._attributes['name'].split('/')[0].split('-')[2])
        self._state = self._status['state']
        self._total_cost = self._status['cost']
        self._n_jobs = self._status['n_jobs']
        self._n_succeeded = self._status['n_succeeded']
        self._n_failed = self._status['n_failed']
        self._n_cancelled = self._status['n_cancelled']
        self._n_completed = self._status['n_completed']

        attempts = await bounded_gather(*[partial(job.attempts) for job in self.jobs], cancel_on_error=True, parallelism=10)
        self._attempts = [att for att in attempts if att is not None]

        child_job_groups = {(await child_jg.attributes()).get('name', ''): child_jg async for child_jg in self.job_group.job_groups()}
        child_job_groups = {name.split('/')[1]: child_jg for name, child_jg in child_job_groups.items()}

        phase_job_group = child_job_groups['phase']

        phase_status = await phase_job_group.status()
        ligate_status = await child_job_groups['ligate'].status()

        self._phase_cost = phase_status['cost']
        self._ligate_cost = ligate_status['cost']

        self._other_costs = max(0.0, self._total_cost - self._phase_cost - self._ligate_cost)

        if not self.is_running:
            self._needs_refresh = False

    @property
    def sample_size(self):
        return int(self._attributes['N'])

    @property
    def percent_completed(self):
        return 100 * (self._n_completed / self._n_jobs)

    @property
    def is_running(self):
        return self._n_completed != self._n_jobs

    @property
    def has_started(self):
        return self._n_completed > 0

    @property
    def state(self):
        return self._state

    @property
    def start_time(self):
        if self._attempts:
            self._start_time = min([parse_timestamp_msecs(attempt.get('start_time'))
                                    for job_attempts in self._attempts
                                    for attempt in job_attempts
                                    if attempt is not None])
        if self._start_time:
            return time_msecs_str(self._start_time)
        return None

    @property
    def end_time(self):
        if not self.is_running:
            end_times = [parse_timestamp_msecs(attempt.get('end_time'))
                         for job_attempts in self._attempts
                         for attempt in job_attempts
                         if attempt is not None and attempt.get('end_time') is not None]
            if end_times:
                self._end_time = max(end_times)
        if self._end_time:
            return time_msecs_str(self._end_time)
        return None

    @property
    def duration(self):
        end_time = self._end_time or time_msecs()
        if self._start_time and end_time:
            return humanize_timedelta_msecs(end_time - self._start_time)
        return None

    def phase_duration_stats(self):
        durations = []
        for job_attempts in self._attempts:
            start_time = None
            end_time = None
            for attempt in job_attempts:
                _start_time = parse_timestamp_msecs(attempt.get('start_time')) or time_msecs()
                start_time = _start_time if start_time is None else min(_start_time, start_time)
                _end_time = parse_timestamp_msecs(attempt.get('end_time')) or time_msecs()
                end_time = _end_time if end_time is None else max(_end_time, end_time)
            durations.append(end_time - start_time)

        durations = [duration / 1000 / 60 for duration in durations]
        if not durations:
            return (None, None, None)

        min_duration = min(durations)
        max_duration = max(durations)
        mean_duration = statistics.mean(durations)
        return (min_duration, max_duration, mean_duration)

    def update_table(self, table):
        min_duration_min, max_duration_min, mean_duration_min = self.phase_duration_stats()

        row = [f"{self._sample_group_id}", f"{self._job_group_id}", f"{self.start_time}", f"{self.end_time}",
               f"{self.duration}", f"{self._state}", f"{self._n_jobs}", f"{self._n_completed}",
               f"{self.percent_completed:.2f}", f"{self._n_succeeded}", f"{self._n_failed}", f"{self._n_cancelled}",
               f"{mean_duration_min:.2f}", f"{min_duration_min:.2f}", f"{max_duration_min:.2f}",
               f"${self._phase_cost:.4f}", f"${self._ligate_cost:.4f}", f"${self._other_costs:.4f}", f"${self._total_cost:.4f}",
               f"${self._total_cost / self.sample_size:.4f}"]

        row = [markup_cell(c, self._state if not self.is_running else 'running', self.has_started) for c in row]

        table.add_row(*row)


async def generate_sample_group_table(b: bc.Batch) -> Table:
    global SAMPLE_GROUPS

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

    if len(SAMPLE_GROUPS) == 0:
        job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]
        job_groups = [jg for name, jg in job_groups if name.startswith('sample-group')]
        SAMPLE_GROUPS = await bounded_gather(*[partial(SampleGroupProgress.initialize, b, jg) for jg in job_groups], parallelism=5)
        SAMPLE_GROUPS.sort(key=lambda jg: jg._sample_group_id)

    await bounded_gather(*[partial(sample_group.refresh) for sample_group in SAMPLE_GROUPS], parallelism=5)

    for sample_group in SAMPLE_GROUPS:
        sample_group.update_table(table)

    return table


def get_contig(contig):
    contig = re.sub('^[^0-9]', '', contig)
    try:
        return int(contig)
    except:
        return contig


async def generate_union_table(b: bc.Batch) -> Table:
    table = Table()

    table.add_column("Contig")
    table.add_column("State")

    job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]

    job_groups = [jg for name, jg in job_groups if name == 'union-sample-groups']
    if not job_groups:
        return table

    contig_job_groups = [((await child_jg.attributes())['contig'], child_jg) async for child_jg in job_groups[0].job_groups()]

    contig_job_groups.sort(key=lambda x: get_contig(x[0]))

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
    table.add_column("Duration")

    job_groups = [((await jg.attributes()).get('name', ''), jg) async for jg in b.job_groups()]

    glimpse_cost = sum([(await jg.status())['cost'] for name, jg in job_groups if name.startswith('sample-group')])
    total_cost = (await b.status())['cost']
    union_cost = total_cost - glimpse_cost

    start_time = parse_timestamp_msecs((await b.status())['time_created'])
    end_time = parse_timestamp_msecs((await b.status())['time_completed']) or time_msecs()
    duration = humanize_timedelta_msecs(end_time - start_time)

    row = [f'${total_cost:.2f}', f'${glimpse_cost:.2f}', f'${union_cost:.2f}', f'{duration}']

    table.add_row(*row)

    return table


async def main(batch_id: int):
    batch_client = await bc.BatchClient.create('dummy')
    try:
        b = await batch_client.get_batch(batch_id)
        with Live(refresh_per_second=30) as live:
            while True:
                sample_group_table = await generate_sample_group_table(b)
                union_sample_groups_table = await generate_union_table(b)
                total_costs = await get_total_costs(b)

                panel_group = Group(
                    Panel(sample_group_table,  title='Sample Groups'),
                    Panel(union_sample_groups_table, title='Union Data By Contig'),
                    Panel(total_costs, title='Batch Summary')
                )

                live.update(panel_group)
                await asyncio.sleep(30)
    finally:
        await batch_client.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--batch-id', type=int, required=True)
    args = p.parse_args()

    asyncio.run(main(args.batch_id))
