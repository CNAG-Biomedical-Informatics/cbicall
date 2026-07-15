# PyPI release TODO

## TestPyPI

- [x] Register a separate account at <https://test.pypi.org/account/register/>.
- [ ] Verify the email address.
- [ ] Configure 2FA with a new authenticator.
- [ ] Store the recovery codes and add a second authentication method if possible.
- [x] Create the `testpypi` GitHub environment without repository secrets.
- [x] Register the TestPyPI trusted publisher for:
  - GitHub owner: `CNAG-Biomedical-Informatics`
  - Repository: `cbicall`
  - Workflow: `publish-testpypi.yml`
  - GitHub environment: `testpypi`
- [x] Run the `Publish to TestPyPI` workflow manually. It builds and validates the
  distributions with:

  ```bash
  python3 -m build
  python3 -m twine check dist/*
  ```

- [x] Confirm that the trusted-publishing job uploads to TestPyPI.

- [x] Test the uploaded wheel without changing the system installation:

  ```bash
  python3 -m pip install \
    --index-url https://test.pypi.org/simple/ \
    --no-deps \
    --target /tmp/cbicall-testpypi \
    "cbicall==1.1.0b1"

  PYTHONPATH=/tmp/cbicall-testpypi python3 -m cbicall --version
  PYTHONPATH=/tmp/cbicall-testpypi python3 -m cbicall validate-registry
  PYTHONPATH=/tmp/cbicall-testpypi python3 -m cbicall validate-resources
  ```

TestPyPI is independent from production PyPI. It has separate credentials, does
not reserve the production package name, and may periodically remove projects.

## Production PyPI

- [ ] Wait for production PyPI account recovery.
- [ ] Restore 2FA and store new recovery codes.
- [ ] Configure the `cbicall` trusted publisher for:
  - GitHub owner: `CNAG-Biomedical-Informatics`
  - Repository: `cbicall`
  - Workflow: `publish-pypi.yml`
  - GitHub environment: `pypi`
- [ ] Create the matching `pypi` environment in the GitHub repository.
- [ ] Select and set the intended package version before creating the release.
- [ ] Do not publish the GitHub release until trusted publishing is configured.

Uploaded PyPI files cannot be replaced, even after deletion. If an upload is
incorrect, yank it and publish a new incremented version such as `.dev2` or
`b2`.
