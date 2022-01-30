from email.mime.nonmultipart import MIMENonMultipart
from email.policy import Policy

class MIMEText(MIMENonMultipart):
    def __init__(self, _text: str, _subtype: str = ..., _charset: str | None = ..., *, policy: Policy | None = ...) -> None: ...
