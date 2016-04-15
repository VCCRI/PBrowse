/**
 * Dalliance configuration objects for all the currently prepared genomes
 */

var genomes = {};
genomes["zv9"] = {
	chr:	'1',
	viewStart:	30000000,
	viewEnd:	30030000,
	cookieKey:	'danio-zv9',
	
	coordSystem: {
	speciesName: 'Danio rerio',
	taxon: 7955,
	auth: 'Zv',
	version: 9
	},
	sources:
	[
	 	{
	 		name: 'Genome',
			uri:	'//www.derkholm.net:9080/das/danRer7comp/',
			desc: 'Zebrafish reference genome build Zv9', 
			tier_type: 'sequence',
			provides_entrypoints: true
		},
		{
			name: 'Genes',
			desc: 'Gene structures from Ensembl 61',
			uri:	'//www.derkholm.net:8080/das/dre_61_9a/',
			collapseSuperGroups: true,
			provides_karyotype: true,
			provides_search: true
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from Ensembl 61',
			uri:	'//www.derkholm.net:8080/das/dre_61_9a/',
			stylesheet_uri: '//www.derkholm.net/dalliance-test/stylesheets/ens-repeats.xml'
		}
	],
	searchEndpoint: new DASSource('//www.derkholm.net:8080/das/dre_61_9a/'),
	karyoEndpoint: new DASSource('//www.derkholm.net:8080/das/dre_61_9a/'),
	fullScreen: true,
	disablePoweredBy:true,
	uiPrefix:"//"+window.location.host+"/",
 
	browserLinks: {
		Ensembl: '//www.ensembl.org/Danio_rerio/Location/View?r=${chr}:${start}-${end}',
		UCSC: '//genome.ucsc.edu/cgi-bin/hgTracks?db=danRer7&position=chr${chr}:${start}-${end}'
	}
};

genomes["generic"] = {
	chr:'1',
	viewStart:50000,
	viewEnd:100000,
	cookieKey:'generic',
	coordSystem: {
		speciesName: 'Generic',
		auth: '',
		version: '',
		ucscName: ''
	},
	sources:[
//      {name:'Genome',
//		twoBitURI:"//"+window.location.host+"/files/653",
//		tier_type: 'sequence'},
//		
//		{name: 'Genes',
//		desc: 'Aspergillus Nidulans Genes',
//		bwgURI: "//"+window.location.host+"/files/652",
//		stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode2.xml',
//		collapseSuperGroups: true,
//		trixURI: "//"+window.location.host+"/files/652.ix"
//		}
    ],
	fullScreen: true,
	disablePoweredBy:true,
	uiPrefix:"//"+window.location.host+"/"
};

genomes["hg37"] = {
	chr:	'22',
	viewStart:	30700000,
	viewEnd:	30900000,
	cookieKey:	'human-grc_h37',
	coordSystem: {
		speciesName: 'Human',
		taxon: 9606,
		auth: 'GRCh',
		version: '37',
		ucscName: 'hg19'
	},
	chains: {
		hg18ToHg19: new Chainset('//www.derkholm.net:8080/das/hg18ToHg19/', 'NCBI36', 'GRCh37',
		{
			speciesName: 'Human',
			taxon: 9606,
			auth: 'NCBI',
			version: 36,
			ucscName: 'hg18'
		})
	},
	sources:	
	[
    	{
    		name:	'Genome',
			twoBitURI:	'//www.biodalliance.org/datasets/hg19.2bit',
			tier_type: 'sequence'
		},
		{	
			name: 'Genes',
			desc: 'Gene structures from GENCODE 19',
			bwgURI: '//www.biodalliance.org/datasets/gencode.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode.xml',
			collapseSuperGroups: true,
			trixURI: '//www.biodalliance.org/datasets/geneIndex.ix'
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from Ensembl 59',
			bwgURI: '//www.biodalliance.org/datasets/repeats.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats.xml'
		}
	],
	uiPrefix:"//"+window.location.host+"/",
	fullScreen: true,
	disablePoweredBy:true,
	browserLinks: {
		Ensembl: '//www.ensembl.org/Homo_sapiens/Location/View?r=${chr}:${start}-${end}',
		UCSC: '//genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr${chr}:${start}-${end}',
		Sequence: '//www.derkholm.net:8080/das/hg19comp/sequence?segment=${chr}:${start},${end}'
	},
	hubs: 
	[
	 	'//ngs.sanger.ac.uk/production/ensembl/regulation/hub.txt', 
	 	'//ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt'
 	],
};

genomes["hg38"] = {
	chr:	'22',
	viewStart:	30300000,
	viewEnd:	30500000,
	cookieKey:	'human-grc_h38',

	coordSystem: {
	speciesName: 'Human',
	taxon: 9606,
	auth: 'GRCh',
	version: '38',
	ucscName: 'hg38'
	},
	chains: {
	hg19ToHg38: new Chainset('//www.derkholm.net:8080/das/hg19ToHg38/', 'GRCh37', 'GRCh38',
		{
			speciesName: 'Human',
			taxon: 9606,
			auth: 'GRCh',
			version: 37,
			ucscName: 'hg19'
		})
	},

	sources:
	[
    	{
    		name:	'Genome',
			twoBitURI:	'//www.biodalliance.org/datasets/hg38.2bit',
			tier_type: 'sequence'
		},
		{
			name: 'GENCODE',
			desc: 'Gene structures from GENCODE 21',
			bwgURI: '//www.biodalliance.org/datasets/GRCh38/gencode.v21.annotation.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode2.xml',
			collapseSuperGroups: true,
			trixURI: '//www.biodalliance.org/datasets/GRCh38/gencode.v21.annotation.ix'
		},
		{
			name: 'GENCODEv19',
			disabled: true,
			desc: 'Gene structures from GENCODE 19',
			bwgURI: '//www.biodalliance.org/datasets/gencode.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode.xml',
			collapseSuperGroups: true,
			pennant: '//genome.ucsc.edu/images/19.jpg',
			trixURI: '//www.biodalliance.org/datasets/gene-index.ix',
			mapping: 'hg19ToHg38'
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from UCSC',
			bwgURI: '//www.biodalliance.org/datasets/GRCh38/repeats.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats2.xml'
		}
	],
	uiPrefix:"//"+window.location.host+"/",
	fullScreen: true,
	disablePoweredBy:true,
	hubs: 
	[
       '//ngs.sanger.ac.uk/production/ensembl/regulation/hub.txt', 
       {url: '//ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt', genome: 'hg19', mapping: 'hg19ToHg38'}
	],
};

genomes["ncbi36"] = {
	chr:	'22',
	viewStart:	30000000,
	viewEnd:	30030000,
	cookieKey:	'ndw-human-ncbi36',

	coordSystem: {
		speciesName: 'Human',
		taxon: 9606,
		auth: 'NCBI',
		version: '36',
		ucscName: 'hg18'
	},

	chains: {
	hg19ToHg18: new Chainset('//www.derkholm.net:8080/das/hg19ToHg18/', 'GRCh37', 'NCBI36',
		{
			speciesName: 'Human',
			taxon: 9606,
			auth: 'GRCh',
			version: 37
		})
	},
	sources:	
	[
	    {
	    	name:'Genome',	
			twoBitURI:'//www.biodalliance.org/datasets/hg18.2bit',
			tier_type:'sequence'},
			{name:'Genes',	
			desc:'Gene structures from Ensembl 54',
			uri:'//www.derkholm.net:8080/das/hsa_54_36p/',	
			collapseSuperGroups:true,
			provides_karyotype:true,
			provides_search:true,
			provides_entrypoints:true,
			maxbins:false
		}, 
		{
			name:	'Repeats',	
			uri:		'//www.derkholm.net:8080/das/hsa_54_36p/',	
			stylesheet_uri:	'//www.biodalliance.org/stylesheets/repeats-L1.xml'
		}
	],
	fullScreen: true,
	disablePoweredBy:true,
	uiPrefix:"//"+window.location.host+"/",
	searchEndpoint: new DASSource('//www.derkholm.net:8080/das/hsa_54_36p/'),
	browserLinks: {
		Ensembl: '//ncbi36.ensembl.org/Homo_sapiens/Location/View?r=${chr}:${start}-${end}',
		UCSC: '//genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr${chr}:${start}-${end}',
		Sequence: '//www.derkholm.net:8080/das/hg18comp/sequence?segment=${chr}:${start},${end}'
	}
};

genomes["m38"] = {
	chr:	'19',
	viewStart:	30000000,
	viewEnd:	30100000,
	cookieKey:	'mouse38',

	coordSystem: {
	speciesName: 'Mouse',
	taxon: 10090,
	auth: 'GRCm',
	version: 38,
	ucscName: 'mm10'
	},
	chains: 
	{
		mm9ToMm10: new Chainset('//www.derkholm.net:8080/das/mm9ToMm10/', 'NCBIM37', 'GRCm38',
		{
			speciesName: 'Mouse',
			taxon: 10090,
			auth: 'NCBIM',
			version: 37,
			ucscName: 'mm9'
		})
	}, 
	sources:	
	[
	 	{
	 		name: 'Genome',
			twoBitURI:	'//www.biodalliance.org/datasets/GRCm38/mm10.2bit',
			desc: 'Mouse reference genome build GRCm38',
			tier_type: 'sequence',
			provides_entrypoints: true
		},
		{
			name: 'Genes',
			desc: 'Gene structures from GENCODE M2',
			bwgURI: '//www.biodalliance.org/datasets/GRCm38/gencodeM2.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode.xml',
			collapseSuperGroups: true,
			trixURI: '//www.biodalliance.org/datasets/GRCm38/gencodeM2.ix'
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from UCSC', 
			bwgURI: '//www.biodalliance.org/datasets/GRCm38/repeats.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats2.xml'
		}
	],
	uiPrefix:"//"+window.location.host+"/",
	fullScreen: true,
	disablePoweredBy:true	
};

genomes["m37"] = {
	chr:	'19',
	viewStart:	30000000,
	viewEnd:	30100000,
	cookieKey:	'mouse37',

	coordSystem: {
	speciesName: 'Mouse',
	taxon: 10090,
	auth: 'GRCm',
	version: 38,
	ucscName: 'mm10'
	},
	chains: {
	mm9ToMm10: new Chainset('//www.derkholm.net:8080/das/mm9ToMm10/', 'NCBIM37', 'GRCm38',
		{
			speciesName: 'Mouse',
			taxon: 10090,
			auth: 'NCBIM',
			version: 37,
			ucscName: 'mm9'
		})
	}, 
	sources:	
	[
	 	{
	 		name: 'Genome',
			twoBitURI:	'//www.biodalliance.org/datasets/GRCm38/mm10.2bit',
			desc: 'Mouse reference genome build GRCm38',
			tier_type: 'sequence',
			provides_entrypoints: true
		},
		{
			name: 'Genes',
			desc: 'Gene structures from GENCODE M2',
			bwgURI: '//www.biodalliance.org/datasets/GRCm38/gencodeM2.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode.xml',
			collapseSuperGroups: true,
			trixURI: '//www.biodalliance.org/datasets/GRCm38/gencodeM2.ix'
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from UCSC', 
			bwgURI: '//www.biodalliance.org/datasets/GRCm38/repeats.bb',
			stylesheet_uri: '//www.biodalliance.org/stylesheets/bb-repeats2.xml'
		}
	],
	uiPrefix:"//"+window.location.host+"/",
	fullScreen: true,
	disablePoweredBy:true
};

genomes["ws220"] = {
	chr:	'I',
	viewStart:	3000000,
	viewEnd:	3030000,
	cookieKey:	'worm_ws220',

	coordSystem: {
	speciesName: 'Caenorhabditis elegans',
	taxon: 6239,
	auth: 'WS',
	version: 220 
	},
	sources:
	[
	 	{
	 		name: 'Genome',
			uri:	'//www.derkholm.net:8080/das/cel_61_220/',
			desc: 'Worm reference genome build WS220',
			tier_type: 'sequence',
			provides_entrypoints: true
		},
		{
			name: 'Genes',
			desc: 'Gene structures from Wormbase WS220',
			uri:	'//www.derkholm.net:8080/das/cel_61_220/',
			collapseSuperGroups: true,
			provides_karyotype: true,
			provides_search: true
		},
		{
			name: 'Repeats',
			desc: 'Repeat annotation from Wormbase + Ensembl',
			uri: '//www.derkholm.net:8080/das/cel_61_220/',
			stylesheet_uri: '//www.derkholm.net/dalliance-test/stylesheets/mouse-repeats.xml'
		}
	],
	searchEndpoint: new DASSource('//www.derkholm.net:8080/das/cel_61_220/'),
	karyoEndpoint: new DASSource('//www.derkholm.net:8080/das/cel_61_220/'),
 
	uiPrefix:"//"+window.location.host+"/",
	fullScreen: true,
	disablePoweredBy:true,
	browserLinks: 
	{
		Ensembl: '//www.ensembl.org/Caenorhabditis_elegans/Location/View?r=${chr}:${start}-${end}'
	}
};
