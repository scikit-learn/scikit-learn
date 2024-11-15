import argparse
import os
import shlex
import subprocess


def execute_command(command):
    command_list = shlex.split(command)
    subprocess.run(command_list, check=True, text=True)


def main():
    comment = os.environ["COMMENT"].splitlines()[0].strip()

    # Extract the command-line arguments from the comment
    prefix = "@scikit-learn-bot update lock-files"
    assert comment.startswith(prefix)
    all_args_list = shlex.split(comment[len(prefix) :])

    # Parse the options for the lock-file script
    parser = argparse.ArgumentParser()
    parser.add_argument("--select-build", default="")
    parser.add_argument("--skip-build", default=None)
    parser.add_argument("--select-tag", default=None)
    args, extra_args_list = parser.parse_known_args(all_args_list)

    # Rebuild the command-line arguments for the lock-file script
    args_string = ""
    if args.select_build != "":
        args_string += f" --select-build {args.select_build}"
    if args.skip_build is not None:
        args_string += f" --skip-build {args.skip_build}"
    if args.select_tag is not None:
        args_string += f" --select-tag {args.select_tag}"

    # Parse extra arguments
    extra_parser = argparse.ArgumentParser()
    extra_parser.add_argument("--commit-marker", default=None)
    extra_args, _ = extra_parser.parse_known_args(extra_args_list)

    marker = ""
    # Additional markers based on the tag
    if args.select_tag == "main-ci":
        marker += "[doc build] "
    elif args.select_tag == "scipy-dev":
        marker += "[scipy-dev] "
    elif args.select_tag == "arm":
        marker += "[cirrus arm] "
    elif len(all_args_list) == 0:
        # No arguments which will update all lock files so add all markers
        marker += "[doc build] [scipy-dev] [cirrus arm] "
    # The additional `--commit-marker` argument
    if extra_args.commit_marker is not None:
        marker += extra_args.commit_marker + " "

    execute_command(
        f"python build_tools/update_environments_and_lock_files.py{args_string}"
    )
    execute_command('git config --global user.name "scikit-learn-bot"')
    execute_command('git config --global user.email "noreply@github.com"')
    execute_command("git add -A")
    # Avoiding commiting the scripts that are downloaded from main
    execute_command("git reset build_tools/shared.sh")
    execute_command("git reset build_tools/update_environments_and_lock_files.py")
    execute_command(
        "git reset build_tools/on_pr_comment_update_environments_and_lock_files.py"
    )
    # Using --allow-empty to handle cases where the lock-file has not changed
    execute_command(f'git commit --allow-empty -m "{marker}Update lock files"')
    execute_command("git push")


if __name__ == "__main__":
    main()
