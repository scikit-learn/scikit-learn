import pytest


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
                        'dateutil/Asia/Tokyo', 'dateutil/US/Pacific'])
def tz(request):
    return request.param
