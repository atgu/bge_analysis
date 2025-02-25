import argparse
import hailtop.batch as hb
from typing import List


def register_account(batch_billing_project: str,
                     regions: List[str],
                     remote_tmpdir: str,
                     email: str,
                     first_name: str,
                     last_name: str):
    backend = hb.ServiceBackend(batch_billing_project, regions=regions, remote_tmpdir=remote_tmpdir)
    b = hb.Batch(backend=backend, name='register-batch-service-account-in-terra')
    j = b.new_bash_job()
    j.image('hailgenetics/python-dill:3.11-slim')
    j.command('apt-get update')
    j.command('apt-get install -y git')
    j.command('pip3 install google-auth urllib3 firecloud')
    j.command(f'git clone https://github.com/broadinstitute/terra-tools.git')
    j.command('cd terra-tools/')
    j.command('ls scripts/register_service_account/')
    j.command(f'python3 scripts/register_service_account/register_service_account_no_keyfile.py -e "{email}" -f "{first_name}" -l "{last_name}"')
    b.run()
    backend.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--email', required=True)
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--regions', required=True)
    parser.add_argument('--remote-tmpdir', required=True)
    parser.add_argument('--first-name', required=True)
    parser.add_argument('--last-name', required=True)

    args = parser.parse_args()

    register_account(args.billing_project,
                     [args.regions],
                     args.remote_tmpdir,
                     args.email,
                     args.first_name,
                     args.last_name)
