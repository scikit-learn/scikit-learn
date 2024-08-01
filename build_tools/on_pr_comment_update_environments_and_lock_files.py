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
    args_string = comment[len(prefix) :]  # including any leading spaces

    # Extract the tag from the args (if any)
    parser = argparse.ArgumentParser()
    parser.add_argument("--select-tag", default=None)
    parser.add_argument("--commit-marker", default=None)
    args, _ = parser.parse_known_args(shlex.split(args_string))

    # Determine the marker to add in the commit message (if any)
    marker = ""
    if args.commit_marker is not None:
        marker += args.commit_marker + " "
    # Additional markers based on the tag
    if args.select_tag == "main-ci":
        marker += "[doc build] "
    elif args.select_tag == "scipy-dev":
        marker += "[scipy-dev] "
    elif args.select_tag == "arm":
        marker += "[cirrus arm] "
    elif args_string.strip() == "":
        # No arguments which will update all lock files so add all markers
        marker += "[doc build] [scipy-dev] [cirrus arm] "

    execute_command(
        f"python build_tools/update_environments_and_lock_files.py{args_string}"
    )
    execute_command('git config --global user.name "scikit-learn-bot"')
    execute_command('git config --global user.email "noreply@github.com"')
    execute_command("git add -A")
    execute_command(f'git commit -m "{marker}Update lock files"')
    execute_command("git push")


if __name__ == "__main__":
    main()
