#!/usr/bin/env python

import os
import subprocess
import shutil

def makedirs(dirname):
	try:
		os.makedirs(dirname)
	except:
		pass

def run(args, env=None):
	process = subprocess.Popen(args, env=env)
	process.wait()

def copydir(src, dst):
	try:
		shutil.rmtree(dst)
	except:
		pass
	shutil.copytree(src, dst)

def copyfile(src, dst):
	shutil.copy2(src, dst)

def compile():
	makedirs('build')
	
	srcDir = os.path.join('src', 'antigen')
	srcFiles = [os.path.join(srcDir, x) for x in os.listdir(srcDir) if x.endswith('.java')]
	args = [
		'javac',
		'-source', '1.7',
		'-target', '1.7',
		'-classpath', os.path.join('lib', '*'),
		'-d', 'build'
	] + srcFiles
	run(args)

def makejar():
	args = [
		'jar',
		'cf',
		os.path.join('dist', 'lib', 'antigen.jar'),
		'-C', 'build',
		'antigen'
	]
	run(args)

if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	compile()
	makedirs('dist')
	copydir('lib', os.path.join('dist', 'lib'))
	makejar()
	copydir('src', os.path.join('dist', 'src'))
	copyfile(os.path.join('scripts', 'run.py'), os.path.join('dist', 'run.py'))
	copydir('example_run', os.path.join('dist', 'example_run'))
	copydir('example_sweep', os.path.join('dist', 'example_sweep'))

