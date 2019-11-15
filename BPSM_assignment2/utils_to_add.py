

def choice_format(input_string):
    pos_responses = ['y', 'yes', 'true']
    neg_responses = ['n', 'no', 'false']
    quit_responses = ['q', 'quit', 'terminate']
    cont_responses = ['', '.']
    input_string = input_string.lower().strip()
    if input_string in pos_responses: 
        return True
    elif input_string in neg_responses:
        return False
    elif input_string in quit_responses: 
        sys.exit(0)
    elif input_string in cont_responses:
        return 'default'
    else: 
        raise Exception('Input string:'+str(input_string)+' was not in a format recognised for parsing')


def check_or_default(value, default, prompt):
    """
    So that the user can specify a seperate value of a var if prompt = True 
    or if prompt = False, uses eiher specified value or default for all functions
    which use that value
    """
    #DO SOMETHING OR IGNORE?
    pass