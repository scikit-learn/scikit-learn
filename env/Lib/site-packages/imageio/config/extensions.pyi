from typing import List, Dict, Optional

class FileExtension:
    extension: str
    priority: List[str]
    name: Optional[str] = None
    description: Optional[str] = None
    external_link: Optional[str] = None
    volume_support: bool

    def __init__(
        self,
        *,
        extension: str,
        priority: List[str],
        name: str = None,
        description: str = None,
        external_link: str = None
    ) -> None: ...
    def reset(self) -> None: ...

extension_list: List[FileExtension]
known_extensions: Dict[str, List[FileExtension]]
video_extensions: List[FileExtension]
