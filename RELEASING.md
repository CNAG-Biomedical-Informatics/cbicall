# Releasing CBIcall

Git tags are the canonical release records. Stable publication does not require
a GitHub Release object.

## Stable release

1. Set the canonical Python version and add its dated `CHANGELOG.md` entry.
2. Run the full release tests.
3. Commit and push `main`.
4. Create and push an annotated tag:

   ```bash
   git tag -a vX.Y.Z -m "Tagging version X.Y.Z" <commit>
   git push origin vX.Y.Z
   ```

5. Confirm that `publish-pypi.yml` succeeds.
6. Run the Docker workflow manually with `vX.Y.Z` as its tag input.

The production workflow accepts only annotated stable tags whose `vX.Y.Z`
value matches the Python package version and `CHANGELOG.md` entry. PyPI Trusted
Publishing must authorize `publish-pypi.yml` with the `pypi` environment.

Never reuse a published version, move a published tag, or delete a tag when
removing a GitHub Release object.

## TestPyPI prerelease

Set a prerelease Python version and manually dispatch `publish-testpypi.yml`
from `main`. TestPyPI publication does not create or require a Git tag.
