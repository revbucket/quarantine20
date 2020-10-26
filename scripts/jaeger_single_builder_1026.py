""" Back on my bullshit 
	This time I'm going to run many iterations for:
	+ All real networks (full sweep) with R=0.5
	+ Vary R for Fb.artist 
	+ Looots of iterations for BA10, 10k 

"""
import sys
sys.path.append('..')
import graph_generators as gg


def jaeger_real_line_rHalf(g_style, g_param, name):
	""" Generates a single jaeger script line """
	base = 'python jaeger_single_real.py ' 
	args = []
	args.append(' -g %s' % g_style )
	args.append(' -gparam %s' % g_param)
	args.append(' -tau 0.5')
	args.append('-gamma 1')
	args.append(' -init 10')
	args.append(' -iter 25')
	args.append(' -name %s' % name)

	script_str = base + ' '.join(args)
	screen = 'screen -dmS %s ' % name 
	return screen + script_str


def build_real_rHalf():
	# Builds 
	base = 'python jaeger_single_real.py' 
	rows = []
	keys = {'gemsec_fb': ['artist', 'athletes', 'company', 'government', 
	                      'new_sites'], 
		    'gemsec_deezer': gg.GEMSEC_DEEZER_ARGS, 
		    'arxiv': gg.ARXIV_COLLAB_ARGS}
	for k, v in keys.items():
		for vel in v:
			rows.append(jaeger_real_line_rHalf(k, vel, '_'.join([k, vel, 'rHalf'])))
	return rows 


def build_fbArtist_rsweep():
	base = 'python jaeger_single_real.py'

	rows = []
	for i in range(-5, 6, 1):
		R = 2 ** i
		args = [' -g gemsec_fb', '-gparam artist', '-gamma 1', '-init 10', 
				'-iter 25']		
		name = 'gemsec_fb_ARTIST_R_%s' % R
		args.append('-tau %s' % R)
		args.append('-name %s' % name)
		rows.append('screen -dmS %s ' % name + base + ' '.join(args))
	return rows 


def build_baseline_bigiter():
	ba10 = ('python jaeger_single_synthetic.py -g ba -gparam 10 -tau 0.5 ' + 
			'-gamma 1.0 -init 10 -iter 100 -name ba10_bigiter')
	fbartist = ('python jaeger_single_real.py -g gemsec_fb -gparam artist ' + 
		        '-tau 0.5 -gamma 1.0 -init 10 -iter 100 -name fbartist_bigiter')
	return ['screen -dmS ' + _ for _ in [ba10, fbartist]]


if __name__ == '__main__':
	rows = build_real_rHalf() + build_fbArtist_rsweep() + build_baseline_bigiter()








