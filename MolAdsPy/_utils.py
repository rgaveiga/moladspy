from functools import wraps


def apply_owner_rules(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        owner = self._belongs_to
        call_handler = kwargs.get("call_handler", True)

        owner._begin_handler(
            self, func, *args, **kwargs
        ) if owner is not None and call_handler else None
        func(self, *args, **kwargs)
        owner._end_handler(
            self, func, *args, **kwargs
        ) if owner is not None and call_handler else None

    return wrapper
