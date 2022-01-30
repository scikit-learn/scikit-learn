from typing import ClassVar, Tuple

class struct_passwd(Tuple[str, str, int, int, str, str, str]):
    pw_name: str
    pw_passwd: str
    pw_uid: int
    pw_gid: int
    pw_gecos: str
    pw_dir: str
    pw_shell: str

    n_fields: ClassVar[int]
    n_sequence_fields: ClassVar[int]
    n_unnamed_fields: ClassVar[int]

def getpwall() -> list[struct_passwd]: ...
def getpwuid(__uid: int) -> struct_passwd: ...
def getpwnam(__name: str) -> struct_passwd: ...
