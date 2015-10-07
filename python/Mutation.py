class Mutation:
    """Store information about a single mutation event.

    This class packs information about a single mutation into 
    a single place for ease of access and modularity. In an 
    effort to save space, individual fields that will appear 
    eventually in Genomediff files are conglomerated into 
    a single fiels (string) to save space, since different 
    types of mutations must print different information during 
    output."""
    def __init__(self, mut_type=None, location=None, length=1, string=None):
        self.mutation_type = mut_type
        self.length = length
        self.string = string
        self.location = location
        return

    """Getters"""
    def get_string(self):
        return self.string
    def get_mutation_type(self):
        return self.mutation_type
    def get_length(self):
        return self.length
    def get_location(self):
        return self.location

    """Setters"""
    def set_string(self, str):
        self.string = str
        return
    def set_length(self, length):
        self.length = length
        return
    def set_type(self, typ):
        self.mutation_type = typ
        return
    def set_loc(self, loc):
        self.location = loc
        return
        
