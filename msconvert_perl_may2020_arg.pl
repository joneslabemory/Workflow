#!/usr/bin/perl
#Author: Karan Uppal
#Date: 05/18/2020

#1) raw files location F:\predictHD
$raw_files_loc=$ARGV[0];

#2) CDF/mzXML file location
$cdf_files_loc=$ARGV[1]; #"G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\_Orbitrap_RAW\\Velos_Raw\\ViLinh\\fridovich_amnioticfluid_serum2015\\serum\\c18pos\\";

#3) Conversion mode: centroid for xcms; profile for apLCMS
$conversion_mode=$ARGV[2]; #"profile"; #options: "centroid" or "profile"

#4)Number of batches
$numbatches=$ARGV[3];

#5)directory structure: hierarchical (batch1, batch2,batch3,...batchn) or regular (all files in one folder)
$dirstructure=$ARGV[4]; #"hierarchical";

#6) select files
$select_files_name=$ARGV[5]; #"G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\_Orbitrap_RAW\\Velos_Raw\\ViLinh\\fridovich_amnioticfluid_serum2015\\serum\\filenames_c18pos.txt"; #"F:\\predictHD\\c18posfilenames.txt";

#7) proteowizard location
$prot_wiz_loc=$ARGV[6]; #"C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.18264.d4d63fbdd\\";
############################No changes needed below this line###########################
#change directory to raw files location
$cmd1="cd ".$raw_files_loc;
system($cmd1);

%filenames = ();

print $select_files_name,"\n";

if(length($select_files_name)>1){
	
open(inf1, "<".$select_files_name);
while($line=<inf1>){
	chomp $line;
	
	$line=~ s/"//g;
	my @columns = split (/\t/,$line);
        my $line    = $columns[0];
	$line=$line.".raw";
	
	$filenames{$line}=1;
}
}
#open new file in the raw files folder to write file names
open(outf, ">Rawfiles.txt");
#chmod 0755,"Rawfiles.txt";

for($i=1;$i<=$numbatches;$i++){
	
	if($dirstructure eq "hierarchical"){
	$raw_files_locsub=$raw_files_loc."\\batch".$i;
	}else{
			$raw_files_locsub=$raw_files_loc;
	}
	
	
	
#get file names of the files in the folder
opendir(DIR, $raw_files_locsub) or die $!;

    while (my $file = readdir(DIR)) {
		
		if ($file =~ m/\.raw/){
			
			if(length($select_files_name)>0){
				if(exists($filenames{$file}))
				{ 
					
					print $file,"\n";
					print outf "$raw_files_locsub\\$file\n";
				}
			}else{
				print $file,"\n";
				print outf "$raw_files_locsub\\$file\n";
			}
		}
    }
    closedir(DIR);
}
    close(outf);
    

if($conversion_mode eq "centroid"){ #\"C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.7398\\msconvert.exe\"
#Execute MSConvert in centroid mode with PeakPicking option set to true
$cmd="'$prot_wiz_loc\msconvert.exe'-f Rawfiles.txt -o '$cdf_files_loc' --64 --zlib --mzXML --filter \"peakPicking true 1-\"";

}else{
	#Execute MSConvert
	$cmd="\"$prot_wiz_loc\msconvert.exe\" -f Rawfiles.txt -o $cdf_files_loc --64 --zlib --mzXML";
	#$cmd="\"C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.18264.d4d63fbdd\\msconvert.exe\" -f Rawfiles.txt -o '$cdf_files_loc' --64 --zlib --mzXML";
	
}
print $cmd;
system($cmd);