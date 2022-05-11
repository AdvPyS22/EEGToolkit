"""
Auxiliary functions to work within the streamlit environment
"""

from EEGData import EEGData, supported_filetypes
import streamlit as st 
from copy import deepcopy
import os 

class Session:
    """
    A class to handle storing values into and getting them from the streamlit session state.
    """
    def __init__(self):

        # a list of keys to include in the 
        # export session dictionary
        self._to_export = set()

    def add( self, key, value, export = True ):
        """
        Add a new item to the session but 
        can also change an existing item.

        Parameters
        ----------
        key : str
            The key of the new item
        value 
            The value to store of the new item
        export : bool
            Wether or not to include the item in the exported session dictioanry.
        """
        st.session_state[ key ] = value

        if export: 
            self._to_export.add( key )
    
    def set( self, key, value, export = True ):
        """
        A synonym of `add`
        """
        self.add( key, value, export )

    def remove( self, key ):
        """
        Remove an item from the session
        """
        try: 
            st.session_state.pop( key )
            self._to_export.remove( key )
        except: 
            pass

    def get( self, key, default = None, rm = False, copy = False ):
        """
        Tries to get an item from the session 
        and returns the default if this fails. 
        If rm = True is set, the item will be 
        removed after getting. 
        If copy = True is set, it will return a deepcopy.
        """
        try: 
            value = st.session_state[ key ]
            if copy: 
                value = deepcopy(value)
        except: 
            value = default

        if rm: 
            self.remove( key )

        return value 

    def setup( self, key, value ):
        """
        Stores an item to the session only if it's not already present.
        """
        if not self.exists( key ):
            self.add( key, value )

    def exists( self, key ):
        """
        Checks if an item is in the session. Returns True if so.
        """
        verdict = key in st.session_state
        return verdict
    
    def hasvalue( self, key ):
        """
        Checks if an item is in the session and is also not None.
        Returns True if both conditions are fulfilled.
        """
        exists = self.exists( key )
        not_none = self.get( key, None ) is not None
        verdict = exists and not_none
        return verdict

    def export( self ):
        """
        Exports a copy of the session state for better readablity wherein non-intelligible entries are removed. 
        """
        to_export = st.session_state[ self._to_export ]
        return to_export



class stEEGData( EEGData ):
    """
    A streamlit compatible version of EEGData
    The primary difference here is that it replaces the original filepath variables
    with filepath.name arguments to work with the string and not the UploadedFile objects.
    """
    def __init__( self, *args, **kwargs ):
        super().__init__(*args, **kwargs)
        
    def _filesuffix(self, filepath):
        """
        Returns the suffix from a filepath
        """
        suffix = os.path.basename( filepath.name )
        suffix = suffix.split(".")[-1]
        return suffix

    def _csv_delimiter(self, filepath):
        """
        Checks if a csv file is , or ; delimited and returns the 
        correct delmiter to use...
        """
        # open the file and read 
        content = deepcopy(filepath).read().decode()

        # check if a semicolon is present
        # if so, we delimit at ; 
        has_semicolon = ";" in content
        delimiter = ";" if has_semicolon else ","

        return delimiter
    
    def _check_sanity(self, signal_path, event_path, sampling_frequency):
        """
        Checks if valid data inputs were provided
        """

        # check if the datafiles conform to suppored filetypes
        fname = os.path.basename(signal_path.name)
        if not any( [ fname.endswith(suffix) for suffix in supported_filetypes ] ):
            suffix = fname.split(".")[-1]
            raise TypeError( f"The signal datafile could not be interpreted ('.{suffix}'), only {supported_filetypes} files are supported!" )

        fname = os.path.basename(event_path.name)
        if not any( [ fname.endswith(suffix) for suffix in supported_filetypes ] ):
            suffix = fname.split(".")[-1]
            raise TypeError( f"The event datafile could not be interpreted ('.{suffix}'), only {supported_filetypes} files are supported!" )
        
        # check if the frequency is a positive number
        if not sampling_frequency > 0:
            raise ValueError( f"Sampling frequency {sampling_frequency} must be a positive value" )