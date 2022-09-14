from ftplib import FTP
ftp = FTP('ftp.ensembl.org')
ftp.login()
ftp.cwd('/pub/release-95/variation/vcf/') 

directories = ftp.nlst()
for directory in directories:
	print 'Entering {}'.format(directory)
	ftp.cwd(directory)
	for filename in  ftp.nlst():
		if filename.find('consequences') != -1:
			print 'Downloading {}'.format(filename)
#                with open(filename, 'wb') as f:
#                    ftp.retrbinary('RETR ' + filename, f.write)
	ftp.cwd('../')
