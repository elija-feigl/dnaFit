import os
import ipywidgets as widgets

class FileBrowser(object):
    """
    https://gist.github.com/DrDub/6efba6e522302e43d055
    """
    def __init__(self):
        self.path = os.getcwd()
        self._update_files()
        
    def _update_files(self):
        self.files = list()
        self.dirs = list()
        if(os.path.isdir(self.path)):
            for f in os.listdir(self.path):
                ff = self.path + "/" + f
                if os.path.isdir(ff):
                    self.dirs.append(f)
                else:
                    self.files.append(f)
        
    def widget(self):
        box = widgets.VBox()
        self._update(box)
        return box
    
    def _update(self, box):

        layout_button = widgets.Layout(width='800px', height='30px', text_align="start", border="1px solid black")

        def on_click(b):
            if b.description == '..':
                self.path = os.path.split(self.path)[0]
            else:
                self.path = self.path + "/" + b.description
            self._update_files()
            self._update(box)
        
        buttons = []
        if self.files:
            button = widgets.Button(description='..', background_color='#d0d0ff',  layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)
        for f in self.dirs:
            button = widgets.Button(description=f, background_color='#d0d0ff',  layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)
        #for f in self.files:
        #    button = widgets.Button(description=f)
        #    button.on_click(on_click)
        #    buttons.append(button)
        box.children = tuple([widgets.HTML("<h3>%s</h3>" % (self.path,))] + buttons)

