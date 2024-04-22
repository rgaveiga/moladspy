from functools import wraps


def apply_owner_rules(obj):
    def outer_wrapper(func):
        @wraps(func)
        def inner_wrapper(*args, **kwargs):
            owner = obj._belongs_to

            owner._begin_handler(obj, func, **kwargs) if owner is not None else None
            func(*args, **kwargs)
            owner._end_handler(obj, func, **kwargs) if owner is not None else None

        return inner_wrapper

    return outer_wrapper
