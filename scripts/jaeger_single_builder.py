""" Look at this bullshit: 
	I'm gonna use this python script to build a bash script that will then run remotely
"""

def jaeger_synth_line(g_style, tau, name, g_params):
	""" Generates a single jaeger script line """
	base = 'python jaeger_single_synthetic.py ' 
	args = []
	args.append(' -g %s' % g_style )
	args.append(' -N 1e4')
	args.append(' -gparam ' + ' '.join(str(_) for _ in g_params))
	args.append(' -tau %s' % tau)
	args.append('-gamma 1')
	args.append(' -init 10')
	args.append(' -iter 25')
	args.append(' -name %s' % name)

	script_str = base + ' '.join(args)
	screen = 'screen -dmS %s ' % name 
	return screen + script_str

def ba_namer(m, tau):
	return 'ba_%s_%s' % (m, tau)

def plc_namer(m, p, tau):
	return 'plc_%s_%s_%s' % (m, p, tau)

def rw_namer(qe,qv,tau):
	return 'rw_%s_%s_%s' % (qe, qv, tau)

def nn_namer(u, k, tau):
	return 'nn_%s_%s_%s' % (u, k, tau)

def ws_namer(k, p, tau):
	return 'ws_%s_%s_%s' % (k, p, tau)


def main():

	BA_PARAMS = [(5,), (10,)] 
	PLC_PARAMS = [(10, .25),(5, .5)] 
	RW_PARAMS = [(.91, .94)] 
	NN_PARAMS = [(.88, 6)]
	WS_PARAMS = [(10, .05)]

	lines = []
	for i in range(-5, 4):
		tau = (2 ** i)
		for ba_param in BA_PARAMS:
			name = ba_namer(*ba_param, tau)
			lines.append(jaeger_synth_line('ba', tau, name, ba_param))
		for plc_param in PLC_PARAMS:
			name = plc_namer(*plc_param, tau)
			lines.append(jaeger_synth_line('plc', tau, name, plc_param))
		for rw_param in RW_PARAMS:
			name = rw_namer(*rw_param, tau)
			lines.append(jaeger_synth_line('rw', tau, name, rw_param))
		for nn_param in NN_PARAMS:
			name = nn_namer(*nn_param, tau)
			lines.append(jaeger_synth_line('nn', tau, name, nn_param))

		for ws_param in WS_PARAMS:
			name = ws_namer(*ws_param, tau)
			lines.append(jaeger_synth_line('ws', tau, name, ws_param))

	with open('synth_multiR_1015.sh', 'w') as f:
		f.write('\n'.join(lines))

	print('\n'.join(lines))
	return lines



if __name__ == '__main__':
	lines = main()