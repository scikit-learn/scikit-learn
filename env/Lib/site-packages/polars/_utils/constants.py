from datetime import date, datetime, timezone

# Integer ranges
I8_MIN = -(2**7)
I16_MIN = -(2**15)
I32_MIN = -(2**31)
I64_MIN = -(2**63)
I128_MIN = -(2**127)
I8_MAX = 2**7 - 1
I16_MAX = 2**15 - 1
I32_MAX = 2**31 - 1
I64_MAX = 2**63 - 1
I128_MAX = 2**127 - 1
U8_MAX = 2**8 - 1
U16_MAX = 2**16 - 1
U32_MAX = 2**32 - 1
U64_MAX = 2**64 - 1

# Temporal
SECONDS_PER_DAY = 86_400
SECONDS_PER_HOUR = 3_600
NS_PER_SECOND = 1_000_000_000
US_PER_SECOND = 1_000_000
MS_PER_SECOND = 1_000

EPOCH_DATE = date(1970, 1, 1)
EPOCH = datetime(1970, 1, 1).replace(tzinfo=None)
EPOCH_UTC = datetime(1970, 1, 1, tzinfo=timezone.utc)
