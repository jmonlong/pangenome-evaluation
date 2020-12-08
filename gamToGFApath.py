import json
import fileinput


for line in fileinput.input():
    linej = json.loads(line)
    pathn = []
    for node in linej['path']['mapping']:
        if 'is_reverse' in node['position'] and node['position']['is_reverse']:
            pathn.append(node['position']['node_id'] + '-')
        else:
            pathn.append(node['position']['node_id'] + '+')
    print('P\t{}\t{}\t{}'.format(linej['name'],
                                 ','.join(pathn),
                                 '*,' * (len(pathn)-2) + '*'))
          
