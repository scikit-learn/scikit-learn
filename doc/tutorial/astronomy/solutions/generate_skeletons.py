"""Generate skeletons from the example code"""
import os

exercise_dir = os.path.dirname(__file__)
if exercise_dir == '':
    exercise_dir = '.'

skeleton_dir = os.path.abspath(os.path.join(exercise_dir,'../skeletons/'))
if not os.path.exists(skeleton_dir):
    os.makedirs(skeleton_dir)

L = os.listdir(exercise_dir)

for f in L:
    if not f.endswith('.py'):
        continue
    
    if f == os.path.basename(__file__):
        continue

    print "parsing %s" % f

    input_file = open(os.path.join(exercise_dir, f))
    output_file = open(os.path.join(skeleton_dir, f), 'w')

    in_exercise_region = False
    
    for line in input_file:
        linestrip = line.strip()
        if linestrip.startswith('#{{{'):
            in_exercise_region = True
            message = linestrip.lstrip('#{{{')
            output_file.write(line.split('#')[0] + '# TODO: %s\n' % message)
        elif in_exercise_region:
            if '#}}}' in line:
                in_exercise_region = False
        else:
            output_file.write(line)

    output_file.close()
