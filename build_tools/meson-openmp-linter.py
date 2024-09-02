# %%
import json
import re
import subprocess
from pathlib import Path


def uses_openmp(targets):
    # TODO is the flag always openmp in all settings maybe try macOS just to
    # make sure (Linux -fopenmp, Windows /openmp should work) ???
    use_openmp = False
    for target in targets["target_sources"]:
        if any("openmp" in arg for arg in target["parameters"]):
            use_openmp = True

    return use_openmp


def get_targets_using_openmp(json_content):
    targets = [target for target in json_content if uses_openmp(target)]
    for target in targets:
        sources = target["target_sources"]
        # TODO is it always in this order (compiler first and linker second)?
        # Let's say probably and revisit in case of issues
        # TODO raise more user-friendly ValueError
        assert "compiler" in sources[0]
        assert "linker" in sources[1]

    return targets


def get_canonical_name_meson(target, build_path):
    # TODO taking [0] I don't expect to have 2 filenames here
    path = Path(target["filename"][0])
    rel_path = path.relative_to(build_path.absolute())
    # Needed on Windows to match git grep output
    rel_path = rel_path.as_posix()
    # OS-specific pattern .cpython on Linux and something like cp312 on Windows
    pattern = r"\.(cpython|cp\d+).+"
    return re.sub(pattern, "", str(rel_path))


def get_canonical_name_git_grep(filename):
    return re.sub(r"\.pyx(\.tp)?", "", filename)


def get_meson_info():
    build_path = Path("build/introspect")
    subprocess.check_call(["meson", "setup", build_path, "--reconfigure"])

    json_out = subprocess.check_output(
        ["meson", "introspect", build_path, "--targets"], text=True
    )
    json_content = json.loads(json_out)
    meson_targets = get_targets_using_openmp(json_content)
    return [get_canonical_name_meson(each, build_path) for each in meson_targets]


def get_git_grep_info():
    git_grep_filenames = (
        subprocess.check_output(
            ["git", "grep", "-lP", "cython.*parallel|_openmp_helpers"], text=True
        )
        .splitlines()
    )
    git_grep_filenames = [f for f in git_grep_filenames if ".pyx" in f]

    return [get_canonical_name_git_grep(each) for each in git_grep_filenames]


def main():
    from_meson = set(get_meson_info())
    from_git_grep = set(get_git_grep_info())

    only_in_git_grep = from_git_grep - from_meson
    only_in_meson = from_meson - from_git_grep

    msg = ""
    if only_in_git_grep:
        only_in_git_grep_msg = "\n".join(
            [f"  {each}" for each in sorted(only_in_git_grep)]
        )
        msg += (
            "Some Cython files use OpenMP,"
            " but their meson.build is missing the openmp_dep dependency,"
            " see below.\n"
            f"{only_in_git_grep_msg}\n\n"
        )

    if only_in_meson:
        only_in_meson_msg = "\n".join([f"  {each}" for each in sorted(only_in_meson)])
        msg += (
            f"Some Cython files do not use OpenMP,"
            " so their meson.build does not need a openmp_dep dependency:\n"
            f"{only_in_meson_msg}\n\n"
        )

    if from_meson != from_git_grep:
        raise ValueError(
            f"Some issues have been found in Meson OpenMP dependencies:\n\n{msg}"
        )


if __name__ == "__main__":
    main()
