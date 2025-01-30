def parse_line_num_spec(spec: str, total: int) -> list[int]:
    """Parse a line number spec (such as "1,2,4-6") and return a list of
    wanted line numbers.
    """
    items = []
    parts = spec.split(',')
    for part in parts:
        try:
            begend = part.strip().split('-')
            if begend == ['', '']:
                raise ValueError
            if len(begend) == 1:
                items.append(int(begend[0]) - 1)
            elif len(begend) == 2:
                start = int(begend[0] or 1)  # left half open (cf. -10)
                end = int(begend[1] or max(start, total))  # right half open (cf. 10-)
                if start > end:  # invalid range (cf. 10-1)
                    raise ValueError
                items.extend(range(start - 1, end))
            else:
                raise ValueError
        except ValueError as exc:
            msg = f'invalid line number spec: {spec!r}'
            raise ValueError(msg) from exc

    return items
