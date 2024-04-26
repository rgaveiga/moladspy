from functools import wraps


def apply_owner_rules(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        owner = self._belongs_to
        call_handler = kwargs.get("call_handler", True)
        method_args=list(args)

        owner._begin_handler(
            self, func, method_args, **kwargs
        ) if owner is not None and call_handler else None
        func(self, *args, **kwargs)
        owner._end_handler(
            self, func, method_args, **kwargs
        ) if owner is not None and call_handler else None

    return wrapper
