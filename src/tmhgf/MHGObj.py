class block:
    def __init__(self, block_string, alignment):
        self.block_string = block_string
        self.alignment = alignment

    def get_acc(self):
        return self.block_string.split('|')[0]

    def get_start(self):
        return int(self.block_string.split('|')[1])

    def get_end(self):
        return int(self.block_string.split('|')[2])

    def get_direction(self):
        return self.block_string.split('|')[3]

    def get_alignment(self):
        return self.alignment
    
    def get_sequence(self):
        return self.alignment.replace("-","").upper()

class mhg:
    def __init__(self, mhg_list):
        self.mhg_list = mhg_list
    
    def add_block(self, block_obj):
        self.mhg_list.append(block_obj)

    def remove_block(self, block_obj):
        self.mhg_list.remove(block_obj)