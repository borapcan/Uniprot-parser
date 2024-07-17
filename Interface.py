import search
import mass_calculator
import isoforms.py
import tryptic cutter

class Interface:
    """ Displays a simple interface when run"""



    def display_menu(self):
        print("""
-------------------------------
---------- Protein searcher------------
-------------------------------\n
Choose options for your protein search:
    1: Protein name
    2: Glycosilation sites
    3: Protein sequence
	4: Protein isoforms 
	5: Tryptic cut of the protein
	5: Protein mass calculator
	6: Protein isoform mass calculator
    7: Protein tryptic cut mass calculator

To see your notes use: show

To remove use: remove <note name>

To exit use: exit
""")
    
	def _note_maker(self,lst):
        if :
            print("Invalid input! Separate words with colons! (':')")
        else:
            try:
                self.nb.make_task(*lst)
            except ValueError:
                print('Ooooopsies, you made a boo-boo!')	

    def run(self):
        while True:
            self.display_menu()
            note_in = input('Write your note: ').split(":")
            if 'exit' in note_in:
                self.nb.quit()
                sys.exit(0)
            elif 'show' in note_in:
                if len(note_in) == 2:
                    self.nb.show(note_in[1])
                    sleep(5)
                else:
                    self.nb.show()
                    sleep(5)
            elif 'remove' in note_in:
                if len(note_in) == 2:
                    print("Attempting remove...")
                    print(note_in[1])
                    sleep(2)
                    self.nb.rem_task(note_in[1])
                else:
                    print("Please specify note to remove!")

            else:
                self._note_maker(note_in)

				BIT BUCKET				
