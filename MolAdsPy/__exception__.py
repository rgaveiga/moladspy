class BasicException(Exception):
    def __init__(self,message=""):
        '''
        __init__(message) -> Exception class constructor.

        Parameters
        ----------
        message : string, optional
            Error message to be displayed. The default is "".

        Returns
        -------
        None.
        '''
        self.message=message
      
    def __str__(self):
        return self.message