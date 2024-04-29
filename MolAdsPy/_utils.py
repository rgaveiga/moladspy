from functools import wraps


def update_owner_attributes(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        owner = self._belongs_to
        
        if owner is not None:
            if func.__name__ == "_update":
                call_update=kwargs.get("call_update", True)
                
                func(self, *args, **kwargs)
                
                owner._update(self, **kwargs) if (call_update is not None) else None
            else:
                raise NotImplementedError("This decorator is available only for methods _update() of atom collections!")

        return wrapper

def apply_owner_handlers(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        owner = self._belongs_to
        call_handler = kwargs.get("call_handler", True)
        func_args=list(args)

        owner._begin_handler(
            self, func, func_args, **kwargs
        ) if (owner is not None and call_handler) else None
        
        func(self, *args, **kwargs)

        owner._end_handler(
            self, func, func_args, **kwargs
        ) if (owner is not None and call_handler) else None

    return wrapper
