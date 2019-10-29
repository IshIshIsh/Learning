ssignment 1


## Utility Functions 
def check_type(input, expected, firm = True):
	"""
	Checks the input against an expected type and either throws an error
	if firm != True will allow continuation and return 'default'
	"""
	if isinstance(input, expected) == True:
		return 'pass'
	else:
		realtype = type(input)
		if firm == True:
			raise ValueError('The input: '+str(input)+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
		else: 
			print('WARNING: The input: '+str(input)+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
			print('WARNING: Attempting to return to default value if applicable')
			return 'default'
