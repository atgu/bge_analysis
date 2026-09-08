import asyncio
from typing import Optional, Dict, cast

import hailtop.batch as hb
from hailtop.batch_client.client import Job as SyncJob
from hailtop.batch_client.aioclient import Job, JobGroup as JobGroupBC
from hailtop.batch import job, JobGroup
from hailtop.utils import async_to_blocking
from hailtop.utils.time import parse_timestamp_msecs, time_msecs

class ImputationJobGroup(hb.batch.JobGroup):
    @staticmethod
    def from_job_group(b: 'ImputationJobSubmitter', jg: JobGroupBC) -> 'ImputationJobGroup':
        imp_jg = ImputationJobGroup(b, None, False)
        imp_jg._async_job_group = jg
        return imp_jg

    def __init__(self,
                 batch: 'ImputationJobSubmitter',
                 parent_job_group: Optional['JobGroup'],
                 new_submission: bool,
                 *,
                 attributes: Optional[Dict[str, str]] = None,
                 cancel_after_n_failures: Optional[int] = None
                 ):
        super().__init__(batch, parent_job_group, attributes=attributes, cancel_after_n_failures=cancel_after_n_failures)
        self._new_submission = new_submission

    def get_job(self, job_id: int) -> hb.Job:
        bc_j = Job.submitted_job(self._batch, job_id)
        dummy_j = hb.Job(self._batch, self, 'foo')
        dummy_j._job_id = job_id
        dummy_j._client_job = SyncJob(bc_j)
        dummy_j._submitted = True
        return dummy_j

    def get_or_create_job_group(self, *, attributes: Optional[Dict[str, str]] = None, cancel_after_n_failures: Optional[int] = None) -> 'ImputationJobGroup':
        b = cast(ImputationJobSubmitter, self._batch)
        if attributes is not None:
            if 'name' in attributes:
                maybe_jg = b.get_job_group_by_name(attributes['name'])
                if maybe_jg is not None:
                    return maybe_jg
        return b._create_job_group(self, attributes=attributes, cancel_after_n_failures=cancel_after_n_failures)

    async def new_python_job(self, name: Optional[str] = None, attributes: Optional[Dict[str, str]] = None) -> 'job.PythonJob':
        return await self._batch._new_python_job(self, name=name, attributes=attributes)

    async def new_bash_job(
        self, name: Optional[str] = None, attributes: Optional[Dict[str, str]] = None, shell: Optional[str] = None
    ) -> 'job.BashJob':
        return await self._batch._new_bash_job(self, name=name, attributes=attributes, shell=shell)


class ImputationJobSubmitter(hb.Batch):
    @staticmethod
    def from_batch_id(batch_id: int, **kwargs) -> 'ImputationJobSubmitter':
        return async_to_blocking(ImputationJobSubmitter._async_from_batch_id(batch_id, **kwargs))

    @staticmethod
    async def _async_from_batch_id(batch_id: int, **kwargs) -> 'ImputationJobSubmitter':
        from hailtop.batch.backend import ServiceBackend  # pylint: disable=import-outside-toplevel

        b = ImputationJobSubmitter(**kwargs)
        assert isinstance(b._backend, ServiceBackend), repr(b._backend)
        b._async_batch = await (await b._backend._batch_client()).get_batch(batch_id)
        b._new_submission = False
        return b

    def __init__(self, **batch_kwargs):
        max_jobs_in_flight = batch_kwargs.pop('max_jobs_in_flight')
        ramp_up_minutes = batch_kwargs.pop('ramp_up_minutes')

        super().__init__(**batch_kwargs)

        self.max_jobs_in_flight = max_jobs_in_flight
        self._remaining = max_jobs_in_flight
        self._start_time: Optional[float] = None
        self.ramp_up_msecs = ramp_up_minutes * 60.0 * 1000

        self._remaining = max_jobs_in_flight
        self._gate = asyncio.Event()
        if self._remaining > 0:
            self._gate.set()

        self._submit_loop: Optional[asyncio.Task] = None

        self._new_submission = True
        self._imputation_job_groups = {}

    def create_job_group(self, *, attributes: Optional[Dict[str, str]] = None, cancel_after_n_failures: Optional[int] = None) -> ImputationJobGroup:
        return self._create_job_group(None, attributes=attributes, cancel_after_n_failures=cancel_after_n_failures)

    def _create_job_group(self,
                          parent_job_group: Optional['JobGroup'],
                          *,
                          attributes: Optional[Dict[str, str]] = None,
                          cancel_after_n_failures: Optional[int] = None) -> ImputationJobGroup:
        from hailtop.batch.backend import ServiceBackend  # pylint: disable=import-outside-toplevel
        assert isinstance(self._backend, ServiceBackend)
        jg = ImputationJobGroup(self, parent_job_group, self._new_submission,
                                attributes=attributes, cancel_after_n_failures=cancel_after_n_failures)
        self._job_groups.append(jg)
        return jg

    async def start(self):
        self._imputation_job_groups = await self._get_existing_job_groups()

        if self._submit_loop is None:
            # Record exactly when we started to calculate elapsed time for the ramp-up
            maybe_time_msecs = None
            if self._async_batch is not None:
                status = await self._async_batch.status()
                maybe_time_msecs = parse_timestamp_msecs(status['time_created'])
            self._start_time = maybe_time_msecs or time_msecs()
            self._submit_loop = asyncio.create_task(self._submit_jobs())

    async def _new_python_job(self, job_group: JobGroup, name: Optional[str] = None, attributes: Optional[Dict[str, str]] = None) -> job.PythonJob:
        # The Safe Event Loop
        while True:
            if self._remaining > 0:
                self._remaining -= 1
                # If we just took the very last slot, shut the gate behind us!
                if self._remaining <= 0:
                    self._gate.clear()
                break  # We got our slot, exit the wait loop
            else:
                # No slots left? Wait for the background loop to open the gate again
                await self._gate.wait()

        return super()._new_python_job(job_group, name, attributes)

    async def new_python_job(self, name: Optional[str] = None,
                             attributes: Optional[Dict[str, str]] = None) -> job.PythonJob:
        return await self._new_python_job(self._root_job_group, name=name, attributes=attributes)

    async def _new_bash_job(
        self, job_group: JobGroup, name: Optional[str] = None, attributes: Optional[Dict[str, str]] = None, shell: Optional[str] = None
    ) -> job.BashJob:
        # The Safe Event Loop
        while True:
            if self._remaining > 0:
                self._remaining -= 1
                # If we just took the very last slot, shut the gate behind us!
                if self._remaining <= 0:
                    self._gate.clear()
                break  # We got our slot, exit the wait loop
            else:
                # No slots left? Wait for the background loop to open the gate again
                await self._gate.wait()

        return super()._new_bash_job(job_group, name, attributes, shell)

    async def new_bash_job(
            self, name: Optional[str] = None, attributes: Optional[Dict[str, str]] = None, shell: Optional[str] = None
    ) -> job.BashJob:
        return await self._new_bash_job(self._root_job_group, name=name, attributes=attributes, shell=shell)

    async def get_n_in_flight_jobs(self) -> int:
        status = await self._async_batch.status()
        n_jobs = status['n_jobs']
        n_completed = status['n_completed']
        print(f'n_jobs {n_jobs} n_completed {n_completed}')
        return n_jobs - n_completed

    async def _submit_jobs(self):
        while True:
            try:
                # 1. Submit pending local jobs
                await self._async_run(wait=False, disable_progress_bar=True, delete_scratch_on_exit=False)

                # 2. Calculate dynamic max capacity
                elapsed = time_msecs() - self._start_time
                buffer_msecs = 15 * 60 * 1000
                elapsed_trimmed = max(0, elapsed - buffer_msecs)
                if self.ramp_up_msecs > 0:
                    fraction = min(1.0, elapsed_trimmed / self.ramp_up_msecs)
                    current_max_jobs = max(5, 5 + int(self.max_jobs_in_flight * fraction))
                else:
                    current_max_jobs = self.max_jobs_in_flight

                # 3. Calculate remaining slots
                n_in_flight_jobs = (await self.get_n_in_flight_jobs()) - 1  # don't include the submit job

                # IMPORTANT: If your math is returning crazy high numbers,
                # hard-enforce that remaining doesn't drop below 0 here.
                new_remaining = max(0, current_max_jobs - n_in_flight_jobs)
                self._remaining = new_remaining

                print(
                    f"Max: {current_max_jobs} | In Flight: {n_in_flight_jobs} | Remaining slots updated to: {self._remaining}")

                # 4. Open or close the gate based on the new count
                if self._remaining > 0 and not self._gate.is_set():
                    self._gate.set()
                elif self._remaining <= 0 and self._gate.is_set():
                    self._gate.clear()

                await asyncio.sleep(15)

            except asyncio.CancelledError:
                # 1. This is triggered when we explicitly cancel the task
                print("Graceful shutdown: Submit loop cancelled.")
                break  # Break out of the while True loop!
            except Exception as e:
                print(f"CRITICAL: Error in submit loop: {repr(e)}")
                await asyncio.sleep(15)

    async def close(self):
        """Cleanly shut down the background submit loop."""
        if self._submit_loop and not self._submit_loop.done():
            self._submit_loop.cancel()
            import contextlib
            with contextlib.suppress(asyncio.CancelledError):
                await self._submit_loop
        print("ImputationJobSubmitter shut down successfully.")

    async def _get_existing_job_groups(self) -> Dict[str, 'ImputationJobGroup']:
        job_groups_bc = [jg_bc async for jg_bc in self._async_batch.job_groups()]
        job_groups = {(await jg_bc.attributes())['name']: ImputationJobGroup.from_job_group(self, jg_bc) for jg_bc in job_groups_bc}
        return job_groups

    def get_job_group_by_name(self, name: str) -> 'ImputationJobGroup':
        return self._imputation_job_groups.get(name)

    def get_or_create_job_group(self, *, attributes: Optional[Dict[str, str]] = None, cancel_after_n_failures: Optional[int] = None) -> 'ImputationJobGroup':
        if attributes is not None:
            name = attributes.get('name')
            maybe_jg = self.get_job_group_by_name(name)
            if maybe_jg is not None:
                return maybe_jg
        return self.create_job_group(attributes=attributes, cancel_after_n_failures=cancel_after_n_failures)
