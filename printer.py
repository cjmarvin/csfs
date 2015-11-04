class Printer():
	
	"""
	Print things to stdout on one line dynamically"""
	def __init__(self, data):
		from sys import stdout
		stdout.write("\r\x1b[K" + data.__str__())
		stdout.flush()