# DATA WITH COVID PARAMETERS 

COOPER_KWARGS = ['china', 'sk', 'india', 'aus', 'usa', 'texas', 'italy']
JO_KWARGS =['sk', 'seoul', 'busan', 'daegu'] 

def cooper_params(kwarg='usa'): # MURICA FIRST
	""" Get params from here:
		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7321055/	
		kwarg is which data we have

	~~~~~~~~~ASSUMES #S CAN CHANGE ~~~~~~~~~~~~~~~~


	RETURNS (tau, gamma)
		where tau    := infection rate 
		      gamma  := recovery rate

	"""
	assert kwarg in COOPER_KWARGS

	cooper_dict = {'china': (0.350, 0.035),
	 			   'sk': (0.4, 0.035), 
	 			   'india': (0.200, 0.040), 
	 			   'aus': (0.190, 0.050), 
	 			   'usa': (0.178, 0.015), 
	 			   'texas': (0.140, 0.048),
	 			   'italy': (0.180, 0.037)}
	return cooper_dict[kwarg]


def jo_params(kwarg=None):
	""" Get params from here: 
		https://www.medrxiv.org/content/10.1101/2020.04.13.20063412v1.full.pdf
		kwarg is which data we have (see asserts) 

		will output the reported averages for whichever kwarg we put in
			
		~~~~~~~~~~ ASSUMES tau(t), gamma(t) ARE TIME-DEPENDENT ~~~~~~~~~~~~~

	RETURNS (tau, gamma)
		where tau    := infection rate 
		      gamma  := recovery rate		
	"""	
	assert kwarg in JO_KWARGS

	jo_dict = {'sk': (0.1656, 0.0253), 
	 		   'seoul': (0.0705, 0.0140),
	 		   'busan': (0.0253, 0.0670), 
	 		   'daegu': (0.0191, 0.0387)}

	return jo_dict[kwarg]


