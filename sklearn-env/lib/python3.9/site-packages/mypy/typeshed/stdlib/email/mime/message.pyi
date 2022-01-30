from email.message import Message
from email.mime.nonmultipart import MIMENonMultipart
from email.policy import Policy

class MIMEMessage(MIMENonMultipart):
    def __init__(self, _msg: Message, _subtype: str = ..., *, policy: Policy | None = ...) -> None: ...
