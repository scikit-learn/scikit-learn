from __future__ import annotations

import abc
import importlib.util
import json
import os
import subprocess
import sys
import zoneinfo
from typing import IO, TYPE_CHECKING, Any, Callable, Literal, Optional, TypedDict, Union

if TYPE_CHECKING:
    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias
    from pathlib import Path

from polars._utils.unstable import issue_unstable_warning

# These typedefs are here to avoid circular import issues, as
# `CredentialProviderFunction` specifies "CredentialProvider"
CredentialProviderFunctionReturn: TypeAlias = tuple[
    dict[str, Optional[str]], Optional[int]
]

CredentialProviderFunction: TypeAlias = Union[
    Callable[[], CredentialProviderFunctionReturn], "CredentialProvider"
]


class AWSAssumeRoleKWArgs(TypedDict):
    """Parameters for [STS.Client.assume_role()](https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/sts/client/assume_role.html#STS.Client.assume_role)."""

    RoleArn: str
    RoleSessionName: str
    PolicyArns: list[dict[str, str]]
    Policy: str
    DurationSeconds: int
    Tags: list[dict[str, str]]
    TransitiveTagKeys: list[str]
    ExternalId: str
    SerialNumber: str
    TokenCode: str
    SourceIdentity: str
    ProvidedContexts: list[dict[str, str]]


class CredentialProvider(abc.ABC):
    """
    Base class for credential providers.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    """

    @abc.abstractmethod
    def __call__(self) -> CredentialProviderFunctionReturn:
        """Fetches the credentials."""


class CredentialProviderAWS(CredentialProvider):
    """
    AWS Credential Provider.

    Using this requires the `boto3` Python package to be installed.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    """

    def __init__(
        self,
        *,
        profile_name: str | None = None,
        assume_role: AWSAssumeRoleKWArgs | None = None,
    ) -> None:
        """
        Initialize a credential provider for AWS.

        Parameters
        ----------
        profile_name : str
            Profile name to use from credentials file.
        assume_role : AWSAssumeRoleKWArgs | None
            Configure a role to assume. These are passed as kwarg parameters to
            [STS.client.assume_role()](https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/sts/client/assume_role.html#STS.Client.assume_role)
        """
        msg = "`CredentialProviderAWS` functionality is considered unstable"
        issue_unstable_warning(msg)

        self._check_module_availability()
        self.profile_name = profile_name
        self.assume_role = assume_role

    def __call__(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials for the configured profile name."""
        import boto3

        session = boto3.Session(profile_name=self.profile_name)

        if self.assume_role is not None:
            return self._finish_assume_role(session)

        creds = session.get_credentials()

        if creds is None:
            msg = "unexpected None value returned from boto3.Session.get_credentials()"
            raise ValueError(msg)

        return {
            "aws_access_key_id": creds.access_key,
            "aws_secret_access_key": creds.secret_key,
            "aws_session_token": creds.token,
        }, None

    def _finish_assume_role(self, session: Any) -> CredentialProviderFunctionReturn:
        client = session.client("sts")

        sts_response = client.assume_role(**self.assume_role)
        creds = sts_response["Credentials"]

        expiry = creds["Expiration"]

        if expiry.tzinfo is None:
            msg = "expiration time in STS response did not contain timezone information"
            raise ValueError(msg)

        return {
            "aws_access_key_id": creds["AccessKeyId"],
            "aws_secret_access_key": creds["SecretAccessKey"],
            "aws_session_token": creds["SessionToken"],
        }, int(expiry.timestamp())

    @classmethod
    def _check_module_availability(cls) -> None:
        if importlib.util.find_spec("boto3") is None:
            msg = "boto3 must be installed to use `CredentialProviderAWS`"
            raise ImportError(msg)


class CredentialProviderAzure(CredentialProvider):
    """
    Azure Credential Provider.

    Using this requires the `azure-identity` Python package to be installed.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    """

    def __init__(
        self,
        *,
        scopes: list[str] | None = None,
        storage_account: str | None = None,
        _verbose: bool = False,
    ) -> None:
        """
        Initialize a credential provider for Microsoft Azure.

        This uses `azure.identity.DefaultAzureCredential()`.

        Parameters
        ----------
        scopes
            Scopes to pass to `get_token`
        storage_account
            If specified, an attempt will be made to retrieve the account keys
            for this account using the Azure CLI. If this is successful, the
            account keys will be used instead of
            `DefaultAzureCredential.get_token()`
        """
        msg = "`CredentialProviderAzure` functionality is considered unstable"
        issue_unstable_warning(msg)

        self._check_module_availability()

        self.account_name = storage_account
        # Done like this to bypass mypy, we don't have stubs for azure.identity
        self.credential = importlib.import_module("azure.identity").__dict__[
            "DefaultAzureCredential"
        ]()
        self.scopes = scopes if scopes is not None else ["https://storage.azure.com/"]
        self._verbose = _verbose

        if self._verbose:
            print(
                (
                    "CredentialProviderAzure "
                    f"{self.account_name = } "
                    f"{self.scopes = } "
                ),
                file=sys.stderr,
            )

    def __call__(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials."""
        if self.account_name is not None:
            try:
                creds = {
                    "account_key": self._get_azure_storage_account_key_az_cli(
                        self.account_name
                    )
                }

                if self._verbose:
                    print(
                        "[CredentialProviderAzure]: retrieved account keys from Azure CLI",
                        file=sys.stderr,
                    )
            except Exception as e:
                if self._verbose:
                    print(
                        f"[CredentialProviderAzure]: failed to retrieve account keys from Azure CLI: {e}",
                        file=sys.stderr,
                    )
            else:
                return creds, None  # type: ignore[return-value]

        token = self.credential.get_token(*self.scopes)

        return {
            "bearer_token": token.token,
        }, token.expires_on

    @classmethod
    def _check_module_availability(cls) -> None:
        if importlib.util.find_spec("azure.identity") is None:
            msg = "azure-identity must be installed to use `CredentialProviderAzure`"
            raise ImportError(msg)

    @staticmethod
    def _extract_adls_uri_storage_account(uri: str) -> str | None:
        # "abfss://{CONTAINER}@{STORAGE_ACCOUNT}.dfs.core.windows.net/"
        #                      ^^^^^^^^^^^^^^^^^
        try:
            return (
                uri.split("://", 1)[1]
                .split("/", 1)[0]
                .split("@", 1)[1]
                .split(".dfs.core.windows.net", 1)[0]
            )

        except IndexError:
            return None

    @staticmethod
    def _get_azure_storage_account_key_az_cli(account_name: str) -> str:
        az_cmd = [
            "az",
            "storage",
            "account",
            "keys",
            "list",
            "--output",
            "json",
            "--account-name",
            account_name,
        ]

        cmd = az_cmd if sys.platform != "win32" else ["cmd", "/C", *az_cmd]

        # [
        #     {
        #         "creationTime": "1970-01-01T00:00:00.000000+00:00",
        #         "keyName": "key1",
        #         "permissions": "FULL",
        #         "value": "..."
        #     },
        #     {
        #         "creationTime": "1970-01-01T00:00:00.000000+00:00",
        #         "keyName": "key2",
        #         "permissions": "FULL",
        #         "value": "..."
        #     }
        # ]

        return json.loads(subprocess.check_output(cmd))[0]["value"]


class CredentialProviderGCP(CredentialProvider):
    """
    GCP Credential Provider.

    Using this requires the `google-auth` Python package to be installed.

    .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.
    """

    def __init__(
        self,
        *,
        scopes: Any | None = None,
        request: Any | None = None,
        quota_project_id: Any | None = None,
        default_scopes: Any | None = None,
    ) -> None:
        """
        Initialize a credential provider for Google Cloud (GCP).

        Parameters
        ----------
        Parameters are passed to `google.auth.default()`
        """
        msg = "`CredentialProviderGCP` functionality is considered unstable"
        issue_unstable_warning(msg)

        self._check_module_availability()

        import google.auth
        import google.auth.credentials

        # CI runs with both `mypy` and `mypy --allow-untyped-calls` depending on
        # Python version. If we add a `type: ignore[no-untyped-call]`, then the
        # check that runs with `--allow-untyped-calls` will complain about an
        # unused "type: ignore" comment. And if we don't add the ignore, then
        # he check that runs `mypy` will complain.
        #
        # So we just bypass it with a __dict__[] (because ruff complains about
        # getattr) :|
        creds, _ = google.auth.__dict__["default"](
            scopes=(
                scopes
                if scopes is not None
                else ["https://www.googleapis.com/auth/cloud-platform"]
            ),
            request=request,
            quota_project_id=quota_project_id,
            default_scopes=default_scopes,
        )
        self.creds = creds

    def __call__(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials."""
        import google.auth.transport.requests

        self.creds.refresh(google.auth.transport.requests.__dict__["Request"]())

        return {"bearer_token": self.creds.token}, (
            int(
                (
                    expiry.replace(tzinfo=zoneinfo.ZoneInfo("UTC"))
                    if expiry.tzinfo is None
                    else expiry
                ).timestamp()
            )
            if (expiry := self.creds.expiry) is not None
            else None
        )

    @classmethod
    def _check_module_availability(cls) -> None:
        if importlib.util.find_spec("google.auth") is None:
            msg = "google-auth must be installed to use `CredentialProviderGCP`"
            raise ImportError(msg)


def _maybe_init_credential_provider(
    credential_provider: CredentialProviderFunction | Literal["auto"] | None,
    source: str
    | Path
    | IO[str]
    | IO[bytes]
    | bytes
    | list[str]
    | list[Path]
    | list[IO[str]]
    | list[IO[bytes]]
    | list[bytes]
    | None,
    storage_options: dict[str, Any] | None,
    caller_name: str,
) -> CredentialProviderFunction | CredentialProvider | None:
    from polars.io.cloud._utils import (
        _first_scan_path,
        _get_path_scheme,
        _is_aws_cloud,
        _is_azure_cloud,
        _is_gcp_cloud,
    )

    if credential_provider is not None:
        msg = f"The `credential_provider` parameter of `{caller_name}` is considered unstable."
        issue_unstable_warning(msg)

    if credential_provider != "auto":
        return credential_provider

    if storage_options is not None:
        return None

    verbose = os.getenv("POLARS_VERBOSE") == "1"

    if (path := _first_scan_path(source)) is None:
        return None

    if (scheme := _get_path_scheme(path)) is None:
        return None

    provider = None

    try:
        provider = (
            CredentialProviderAWS()
            if _is_aws_cloud(scheme)
            else CredentialProviderAzure(
                storage_account=(
                    CredentialProviderAzure._extract_adls_uri_storage_account(str(path))
                ),
                _verbose=verbose,
            )
            if _is_azure_cloud(scheme)
            else CredentialProviderGCP()
            if _is_gcp_cloud(scheme)
            else None
        )
    except ImportError as e:
        if verbose:
            msg = f"Unable to auto-select credential provider: {e}"
            print(msg, file=sys.stderr)

    if provider is not None and verbose:
        msg = f"Auto-selected credential provider: {type(provider).__name__}"
        print(msg, file=sys.stderr)

    return provider
