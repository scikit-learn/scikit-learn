from __future__ import annotations

import abc
import importlib.util
import json
import os
import subprocess
import sys
import zoneinfo
from datetime import datetime
from functools import partial
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Optional,
    TypedDict,
    Union,
)

import polars._utils.logging
from polars._utils.logging import eprint, verbose
from polars.io.cloud._utils import NoPickleOption

if TYPE_CHECKING:
    from polars._dependencies import boto3

    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

from polars._utils.unstable import issue_unstable_warning

# These typedefs are here to avoid circular import issues, as
# `CredentialProviderFunction` specifies "CredentialProvider"
CredentialProviderFunctionReturn: TypeAlias = tuple[dict[str, str], Optional[int]]

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


class CachingCredentialProvider(CredentialProvider, abc.ABC):
    """
    Base class for credential providers that has built-in caching.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """

    def __init__(self) -> None:
        self._cached_credentials: NoPickleOption[CredentialProviderFunctionReturn] = (
            NoPickleOption()
        )
        self._has_logged_use_cache = False

    def __call__(self) -> CredentialProviderFunctionReturn:
        if os.getenv("POLARS_DISABLE_PYTHON_CREDENTIAL_CACHING") == "1":
            self._cached_credentials.set(None)

            return self.retrieve_credentials_impl()

        credentials = self._cached_credentials.get()

        if credentials is None or (
            (expiry := credentials[1]) is not None
            and expiry <= int(datetime.now().timestamp())
        ):
            credentials = self.retrieve_credentials_impl()
            self._cached_credentials.set(credentials)
            self._has_logged_use_cache = False

        elif verbose() and not self._has_logged_use_cache:
            expiry = credentials[1]
            eprint(
                f"[{CachingCredentialProvider.__repr__(self)}]: "
                f"Using cached credentials ({expiry = })"
            )
            self._has_logged_use_cache = True

        creds, expiry = credentials

        return {**creds}, expiry

    @abc.abstractmethod
    def retrieve_credentials_impl(self) -> CredentialProviderFunctionReturn: ...

    def clear_cached_credentials(self) -> None:
        self._cached_credentials.set(None)

    def __repr__(self) -> str:
        return f"CachingCredentialProvider[{type(self).__name__} @ {hex(id(self))}]"


class CachedCredentialProvider(CachingCredentialProvider):
    """
    Wrapper that adds caching on top of a credential provider.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """

    def __init__(
        self, provider: CredentialProvider | CredentialProviderFunction
    ) -> None:
        self._provider = provider

        super().__init__()

    def retrieve_credentials_impl(self) -> CredentialProviderFunctionReturn:
        return self._provider()

    def __repr__(self) -> str:
        return f"CachedCredentialProvider[{self._provider!r}]"


class CredentialProviderAWS(CachingCredentialProvider):
    """
    AWS Credential Provider.

    Using this requires the `boto3` Python package to be installed.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """

    def __init__(  # noqa: D417 (TODO)
        self,
        *,
        profile_name: str | None = None,
        region_name: str | None = None,
        assume_role: AWSAssumeRoleKWArgs | None = None,
        _auto_init_unhandled_key: str | None = None,
        _storage_options_has_endpoint_url: bool = False,
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

        self._ensure_module_availability()

        self.profile_name = profile_name
        self.region_name = region_name
        self.assume_role = assume_role
        self._auto_init_unhandled_key = _auto_init_unhandled_key
        self._storage_options_has_endpoint_url = _storage_options_has_endpoint_url

        super().__init__()

    def retrieve_credentials_impl(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials for the configured profile name."""
        assert not self._auto_init_unhandled_key

        session = self._session()

        if self.assume_role is not None:
            return self._finish_assume_role(session)

        creds = session.get_credentials()

        if creds is None:
            msg = "did not receive any credentials from boto3.Session.get_credentials()"
            raise self.EmptyCredentialError(msg)

        expiry = (
            int(expiry.timestamp())
            if isinstance(expiry := getattr(creds, "_expiry_time", None), datetime)
            else None
        )

        return {
            "aws_access_key_id": creds.access_key,
            "aws_secret_access_key": creds.secret_key,
            **({"aws_session_token": creds.token} if creds.token is not None else {}),
        }, expiry

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

    # Called from Rust, mainly for AWS endpoint_url
    def _storage_update_options(self) -> dict[str, str]:
        if self._storage_options_has_endpoint_url:
            return {}

        try:
            config = self._session()._session.get_scoped_config()
        except ImportError:
            return {}

        if endpoint_url := config.get("endpoint_url"):
            if verbose():
                eprint(f"[CredentialProviderAWS]: Loaded endpoint_url: {endpoint_url}")

            return {"endpoint_url": endpoint_url}

        return {}

    # Called from Rust
    def _can_use_as_provider(self) -> bool:
        if self._auto_init_unhandled_key:
            if verbose():
                eprint(
                    "[CredentialProviderAWS]: Will not be used as a provider: "
                    f"unhandled key in storage_options: '{self._auto_init_unhandled_key}'"
                )

            return False

        try:
            self()

        except ImportError as e:
            if self.profile_name:
                msg = (
                    "cannot load requested aws_profile "
                    f"'{self.profile_name}': {type(e).__name__}: {e}"
                )
                raise polars.exceptions.ComputeError(msg) from e

            return False

        except self.EmptyCredentialError:
            if verbose():
                eprint("[CredentialProviderAWS]: Did not find any credentials")

            return False

        return True

    def _session(self) -> boto3.Session:
        # Note: boto3 automatically sources the AWS_PROFILE env var
        import boto3

        return boto3.Session(
            profile_name=self.profile_name,
            region_name=self.region_name,
        )

    @classmethod
    def _ensure_module_availability(cls) -> None:
        if importlib.util.find_spec("boto3") is None:
            msg = "boto3 must be installed to use `CredentialProviderAWS`"
            raise ImportError(msg)

    class EmptyCredentialError(Exception):
        """
        Raised when boto3 returns empty credentials.

        This generally indicates that no credentials could be found in the
        environment.
        """


class CredentialProviderAzure(CachingCredentialProvider):
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
        tenant_id: str | None = None,
        credential: Any | None = None,
        _storage_account: str | None = None,
    ) -> None:
        """
        Initialize a credential provider for Microsoft Azure.

        By default, this uses `azure.identity.DefaultAzureCredential()`.

        Parameters
        ----------
        scopes
            Scopes to pass to `get_token`
        tenant_id
            Azure tenant ID.
        credential
            Optionally pass an instantiated Azure credential class to use (e.g.
            `azure.identity.DefaultAzureCredential`). The credential class must
            have a `get_token()` method.
        """
        msg = "`CredentialProviderAzure` functionality is considered unstable"
        issue_unstable_warning(msg)

        self.account_name = _storage_account
        self.scopes = (
            scopes if scopes is not None else ["https://storage.azure.com/.default"]
        )
        self.tenant_id = tenant_id
        self.credential = credential

        if credential is not None:
            # If the user passes a credential class, we just need to ensure it
            # has a `get_token()` method.
            if not hasattr(credential, "get_token"):
                msg = (
                    f"the provided `credential` object {credential!r} does "
                    "not have a `get_token()` method."
                )
                raise ValueError(msg)

        # We don't need the module if we are permitted and able to retrieve the
        # account key from the Azure CLI.
        elif self._try_get_azure_storage_account_credential_if_permitted() is None:
            self._ensure_module_availability()

        if verbose():
            eprint(
                "[CredentialProviderAzure]: "
                f"{self.account_name = } "
                f"{self.tenant_id = } "
                f"{self.scopes = } "
            )

        super().__init__()

    def retrieve_credentials_impl(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials."""
        if (
            v := self._try_get_azure_storage_account_credential_if_permitted()
        ) is not None:
            return v

        import azure.identity

        credential = self.credential or azure.identity.DefaultAzureCredential()
        token = credential.get_token(*self.scopes, tenant_id=self.tenant_id)

        return {
            "bearer_token": token.token,
        }, token.expires_on

    def _try_get_azure_storage_account_credential_if_permitted(
        self,
    ) -> CredentialProviderFunctionReturn | None:
        POLARS_AUTO_USE_AZURE_STORAGE_ACCOUNT_KEY = os.getenv(
            "POLARS_AUTO_USE_AZURE_STORAGE_ACCOUNT_KEY"
        )

        verbose = polars._utils.logging.verbose()

        if verbose:
            eprint(
                "[CredentialProviderAzure]: "
                f"{self.account_name = } "
                f"{POLARS_AUTO_USE_AZURE_STORAGE_ACCOUNT_KEY = }"
            )

        if (
            self.account_name is not None
            and POLARS_AUTO_USE_AZURE_STORAGE_ACCOUNT_KEY == "1"
        ):
            try:
                creds = {
                    "account_key": self._get_azure_storage_account_key_az_cli(
                        self.account_name
                    )
                }

                if verbose:
                    eprint(
                        "[CredentialProviderAzure]: Retrieved account key from Azure CLI"
                    )
            except Exception as e:
                if verbose:
                    eprint(
                        f"[CredentialProviderAzure]: Could not retrieve account key from Azure CLI: {e}"
                    )
            else:
                return creds, None

        return None

    @classmethod
    def _ensure_module_availability(cls) -> None:
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

    @classmethod
    def _get_azure_storage_account_key_az_cli(cls, account_name: str) -> str:
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

        return json.loads(
            cls._azcli(
                "storage",
                "account",
                "keys",
                "list",
                "--output",
                "json",
                "--account-name",
                account_name,
            )
        )[0]["value"]

    @classmethod
    def _azcli_version(cls) -> str | None:
        try:
            return json.loads(cls._azcli("version"))["azure-cli"]
        except Exception:
            return None

    @staticmethod
    def _azcli(*args: str) -> bytes:
        return subprocess.check_output(
            ["az", *args] if sys.platform != "win32" else ["cmd", "/C", "az", *args]
        )


class CredentialProviderGCP(CachingCredentialProvider):
    """
    GCP Credential Provider.

    Using this requires the `google-auth` Python package to be installed.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.
    """

    def __init__(  # noqa: D417 (TODO)
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

        self._ensure_module_availability()

        import google.auth

        self._init_creds = partial(
            google.auth.default,
            scopes=(
                scopes
                if scopes is not None
                else ["https://www.googleapis.com/auth/cloud-platform"]
            ),
            request=request,
            quota_project_id=quota_project_id,
            default_scopes=default_scopes,
        )

        super().__init__()

    def retrieve_credentials_impl(self) -> CredentialProviderFunctionReturn:
        """Fetch the credentials."""
        import google.auth.transport.requests

        creds, _project_id = self._init_creds()
        creds.refresh(google.auth.transport.requests.Request())  # type: ignore[no-untyped-call, unused-ignore]

        return {"bearer_token": creds.token}, (
            int(
                (
                    expiry.replace(tzinfo=zoneinfo.ZoneInfo("UTC"))
                    if expiry.tzinfo is None
                    else expiry
                ).timestamp()
            )
            if (expiry := creds.expiry) is not None
            else None
        )

    @classmethod
    def _ensure_module_availability(cls) -> None:
        if importlib.util.find_spec("google.auth") is None:
            msg = "google-auth must be installed to use `CredentialProviderGCP`"
            raise ImportError(msg)


class UserProvidedGCPToken(CredentialProvider):
    """User-provided GCP token in storage_options."""

    def __init__(self, token: str) -> None:
        self.token = token

    def __call__(self) -> CredentialProviderFunctionReturn:
        return {"bearer_token": self.token}, None


def _get_credentials_from_provider_expiry_aware(
    credential_provider: CredentialProviderFunction,
) -> dict[str, str] | None:
    if (
        isinstance(credential_provider, CredentialProviderAWS)
        and not credential_provider._can_use_as_provider()
    ):
        return None

    creds, opt_expiry = credential_provider()

    if (
        opt_expiry is not None
        and (expires_in := opt_expiry - int(datetime.now().timestamp())) < 7
    ):
        from time import sleep

        if verbose():
            eprint(f"waiting for {expires_in} seconds for refreshed credentials")

        sleep(1 + expires_in)
        creds, _ = credential_provider()

    # Loads the endpoint_url
    if isinstance(credential_provider, CredentialProviderAWS) and (
        v := credential_provider._storage_update_options()
    ):
        creds = {**creds, **v}

    return creds
