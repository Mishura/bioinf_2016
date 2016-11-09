__author__ = "Mikhail R"
import MySQLdb as mdb
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import argparse

locations = ((1, 'CDS','CDS'),
(2, 'Repeat region','repeat_region'),
(4, 'Ribosomal RNA','rRNA'),
(5, 'Transfer RNA','tRNA'),
(6, 'Miscellaneous RNA','misc_RNA'),
(7, 'Non-coding RNA','ncRNA'),
(8, 'Mobile element','mobile_element'),
(9, 'STS','STS'))

def get_id_location(type):
	for j in locations:
		if type == j[2]:
			return j[0]

def process_file(path_to_genbank):


#####
#	CONNECT TO THE DATABASE
#####

	try:
		cnx = mdb.connect('localhost', 'root', 'mysqlAdmin', 'bioinf_2016');  
		cursor = cnx.cursor()
	except:
		print 'Check correctness of connection to mysql.'

#####
#	Check, whether were gene_list table created ot not.
#####

	try:
		cursor.execute("SELECT * FROM gene_list")
		cursor.fetchone()
	except:
		cursor.execute("""
CREATE TABLE IF NOT EXISTS `gene_list` (
  `gene_id` int(11) NOT NULL AUTO_INCREMENT,
  `gene_name` varchar(45) DEFAULT NULL,
  `locus_tag` varchar(45) DEFAULT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` int(2) NOT NULL,
  `feature_id` int(11) NOT NULL,
  PRIMARY KEY (`gene_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;""")
		cnx.commit()
		print 'It seems that you didn`t create table gene_list in database bioinf_2016. Script did it for your.'


#####
#	TRUNCATE gene_list;
#####

	cursor.execute("TRUNCATE gene_list;");
	cnx.commit()


#####
#	PARSE genbank file and CREATE a query.
#####
	sql = "INSERT INTO `gene_list`(`gene_name`, `locus_tag`, `start`, `end`, `strand`, `feature_id`) VALUES "
	try:
		gb = SeqIO.parse(open(path_to_genbank,"r"), "genbank")
	except:
		print 'THIS FILE DOES NOT EXIST: %s' % path_to_genbank
		return


	for gb_record in gb:
		for ind, feature in enumerate(gb_record.features):
			feature_id = get_id_location(feature.type)
			if feature_id:
				end = feature.location.nofuzzy_end
				start = feature.location.nofuzzy_start
				strand = feature.strand
				try:
					locus_tag = feature.qualifiers.get('locus_tag')[0]
				except:
					locus_tag = None
				try:
					gene_name = feature.qualifiers.get('gene')[0]
				except:
					gene_name = None
				if gene_name is None and locus_tag is not None:
					sql += "(%s,'%s',%s,%s,%s,%s)," % ("NULL", locus_tag, start, end, strand, feature_id)
				elif gene_name is not None and locus_tag is None:
					sql += "('%s',%s,%s,%s,%s,%s)," % (gene_name, "NULL", start, end, strand, feature_id)
				elif gene_name is not None and locus_tag is not None:
					sql += "('%s','%s',%s,%s,%s,%s)," % (gene_name, locus_tag, start, end, strand, feature_id)
				elif gene_name is None and locus_tag is None:
					sql += "(%s,%s,%s,%s,%s,%s)," % ("NULL", "NULL", start, end, strand, feature_id)
				print "INSERT INTO `gene_list`(`gene_name`, `locus_tag`, `start`, `end`, `strand`, `feature_id`) VALUES (%s,%s,%s,%s,%s,%s);" % (gene_name, locus_tag, start, end, strand, feature_id)
	cursor.execute(sql[:-1])
	cnx.commit()
	print 'GREAT! LET`S LOOK AT THE DATABASE bioinf_2016! PLEASE, CONNECT TO MySQL database.'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Parse GenBank file and upload it to the database bioinf_2016 to the table gene_list. If it wasn`t created, it won`t work. ')
	parser.add_argument('path', metavar='<path_to_genbank_file>', type=str,
	                   help='path to GenBank file')
	args = parser.parse_args()
	path_to_genbank = args.path
	print "PARSING FILE", path_to_genbank
	process_file(path_to_genbank)


			


		
