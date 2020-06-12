"""Generate dependency table for docs."""
from io import StringIO
from pathlib import Path

from sklearn._build_utils.dependencies import dependent_packages

REPO_FOLDER = Path(__file__).parent.parent.resolve()

# get length of header
package_header_len = max(len(package) for package in dependent_packages) + 4
version_header_len = len('Minimum Version') + 4
tags_header_len = max(len(tags) for _, tags in dependent_packages.values()) + 4

output = StringIO()
output.write(' '.join(['=' * package_header_len,
                       '=' * version_header_len,
                       '=' * tags_header_len]))
output.write('\n')
dependency_title = "Dependency"
version_title = "Minimum Version"
tags_title = "Purpose"

output.write(f'{dependency_title:<{package_header_len}} '
             f'{version_title:<{version_header_len}} '
             f'{tags_title:<{tags_header_len}}\n')

output.write(' '.join(['=' * package_header_len,
                       '=' * version_header_len,
                       '=' * tags_header_len]))
output.write('\n')

for package, (version, tags) in dependent_packages.items():
    output.write(f'{package:<{package_header_len}} '
                 f'{version:<{version_header_len}} '
                 f'{tags:<{tags_header_len}}\n')

output.write(' '.join(['=' * package_header_len,
                       '=' * version_header_len,
                       '=' * tags_header_len]))
output.write('\n')

output = output.getvalue()
print(output)
with (REPO_FOLDER / "doc" / "dependency.rst").open('w') as f:
    f.write(output)
