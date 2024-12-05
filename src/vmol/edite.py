class EditionTasks:
    """ Function to remove, transport or doing what you want using buttons
    instead of typing functions"""
    def __init__(self, master_object):
        self.vmol = master_object
        self.buttons = []
        
        # Buttons
        self.remo = None
        self.move = None
        self.coox = None
        self.cooy = None
        self.cooy = None
    
    def empty(self):
        return
        
    def clean_start(self):
        self.remo.delete()
        self.move.delete()
        self.remo = None
        self.move = None

        return


    def start(self):
        self.remo = self.vmol.button(text='Delete',
                                     bind=self.delete)

        self.move = self.vmol.button(text='Translate',
                                     bind=self.translate)

    
    def delete(self):
        sel = self.vmol.selected
        if sel is None:
            self.vmol.caption.text = 'Select an object to edit.'
        else:
            self
remove