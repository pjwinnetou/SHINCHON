#!/usr/bin/perl -w

use Cwd;

#$package = "IPGlasma_dAu";
$package = "IPGlasma_ellipse";
$maindir = getcwd();

$groupnum = 67;

open(FILE,"</gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/04.IP-Glasma/source/ipglasma-00/input");
@init_pars = <FILE>;
close(FILE);

$par_length = "15";
$par_size = "256";
$par_sizeOutput = "256";
$par_m = "0.2";
#30 fm
#$par_g2mu = "0.025"; #Qs=1.1 GeV, for 2048
#$par_g2mu = "0.10"; #Qs=1.1 GeV, for 1024
#$par_g2mu = "0.40"; #Qs=1.1 GeV, for 512 
#$par_g2mu = "1.60"; #Qs=1.1 GeV, for 256 
#15 fm
#$par_g2mu = "0.006"; #Qs=1.1 GeV, for 2048
#$par_g2mu = "0.025"; #Qs=1.1 GeV, for 1024
#$par_g2mu = "0.10"; #Qs=1.1 GeV, for 512 
$par_g2mu = "0.40"; #Qs=1.1 GeV, for 256 
$par_useNucleus = "0";
$par_proj = "d";
$par_targ = "Au";
$par_bmin = "0.";
$par_bmax = "10.";
$par_nproj = 2;
$par_ntarg = 197;
$par_quark = 3;
$par_useTimeForSeed = 0;
$par_SigmaNN = "42.";
$par_inverseQsForMaxTime = 0;
$par_useFluctuatingx = 1;
$par_maxtime = 1.0; #default = 0.002

$rundir = "${maindir}/running_grp${groupnum}";
mkdir $rundir;

for ($irun=0; $irun<1000; $irun++){

#	sleep 2;

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
	print FILE "Requirements = CPU_Speed>=2\n";
	print FILE "Rank = CPU_Speed\n";
	print FILE "Getenv = true\n";
	print FILE "Priority = +20\n";
	print FILE "request_memory = 2G\n";
	print FILE "Executable = jobscript\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
	print FILE "Notify_user = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment = \"phenix\"\n";
#	print FILE "+Job_Type = \"cas\"\n";
	print FILE "Queue\n";
	close(FILE);

#	$seednum = int(rand(100000000));

	open(FILE, ">input");
	foreach $init_par (@init_pars){
		if ($init_par=~/^size /){
			$init_par = "size ${par_size}\n";
		}elsif ($init_par=~/^m /){
			$init_par = "m ${par_m}\n";
		}elsif ($init_par=~/^g2mu /){
			$init_par = "g2mu ${par_g2mu}\n";
		}elsif ($init_par=~/^useNucleus /){
			$init_par = "useNucleus ${par_useNucleus}\n";
		}elsif ($init_par=~/^sizeOutput /){
			$init_par = "sizeOutput ${par_sizeOutput}\n";
		}elsif ($init_par=~/^seed /){
			$init_par = "seed ${irun}\n";
		}elsif ($init_par=~/^Projectile /){
			$init_par = "Projectile ${par_proj}\n";
		}elsif ($init_par=~/^Target /){
			$init_par = "Target ${par_targ}\n";
		}elsif ($init_par=~/^bmin /){
			$init_par = "bmin ${par_bmin}\n";
		}elsif ($init_par=~/^bmax /){
			$init_par = "bmax ${par_bmax}\n";
		}elsif ($init_par=~/^L /){
			$init_par = "L ${par_length}\n";
		}elsif ($init_par=~/^LOutput /){
			$init_par = "LOutput ${par_length}\n";
		}elsif ($init_par=~/^useConstituentQuarkProton /){
			$init_par = "useConstituentQuarkProton ${par_quark}\n";
		}elsif ($init_par=~/^useTimeForSeed /){
			$init_par = "useTimeForSeed ${par_useTimeForSeed}\n";
		}elsif ($init_par=~/^useFluctuatingx /){
			$init_par = "useFluctuatingx ${par_useFluctuatingx}\n";
		}elsif ($init_par=~/^SigmaNN /){
			$init_par = "SigmaNN ${par_SigmaNN}\n";
		}elsif ($init_par=~/^inverseQsForMaxTime /){
			$init_par = "inverseQsForMaxTime ${par_inverseQsForMaxTime}\n";
		}elsif ($init_par=~/^maxtime /){
			$init_par = "maxtime ${par_maxtime}\n";
		}
		print FILE $init_par;
	}
	close(FILE);

#	$grpdir = sprintf("%s/%s_grp%03d",$maindir,$package,$groupnum);

#	$randnum = int(rand(300));

	open(FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/sphenix/setup.csh\n\n\n";

#	print FILE "mkdir -p $grpdir\n\n";

	$condor_wrk_dir = "${package}_grp${groupnum}_run${irun}";
	print FILE "mkdir -p \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";
	print FILE "cd \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";
	print FILE "rm -rf *\n";
#
#
	$text_out_dir = "/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/shlim/simulation/IP_Glasma/${package}_grp${groupnum}/run${irun}";
	print FILE "mkdir -p ${text_out_dir}\n";

#	print FILE "mkdir -p run\n";
#	print FILE "cd run\n";

	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/04.IP-Glasma/run/he3_plaintext.dat .\n";
	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/04.IP-Glasma/source/ipglasma-00/* .\n";
	print FILE "make clean\n";
	print FILE "make\n";
	print FILE "rm -rf input\n";
	print FILE "cp -av ${wrkdir}/input .\n";
	print FILE "./ipglasma input\n";

	print FILE "ls -t epsilon* | tac > file.lst\n";
	print FILE "cp -av ${maindir}/draw_def2.C .\n";
	print FILE "root -l -b -q 'draw_def2.C+g($par_length,$par_size,$irun,$par_nproj,$par_ntarg,$par_quark)'\n";

#	$outputfile = sprintf("%s/Epsilon_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$rundir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	$outputfile = sprintf("%s/Epsilon_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$text_out_dir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	print FILE "mv outfile.root ${outputfile}\n"; 

	print FILE "ls -t Tmunu* | tac > file.lst\n";
	print FILE "cp -av ${maindir}/draw_Tmunu.C .\n";
	print FILE "root -l -b -q 'draw_Tmunu.C+g($par_length,$par_size,$irun,$par_nproj,$par_ntarg,$par_quark)'\n";

#	$outputfile = sprintf("%s/Tmunu_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$rundir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	$outputfile = sprintf("%s/Tmunu_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$text_out_dir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	print FILE "mv outfile.root ${outputfile}\n"; 

	print FILE "ls -t NpartdNdy-t* | tac > file.lst\n";
	print FILE "cp -av ${maindir}/draw_dNdy.C .\n";
	print FILE "root -l -b -q 'draw_dNdy.C+g($par_length,$par_size,$irun,$par_nproj,$par_ntarg,$par_quark)'\n";

#	$outputfile = sprintf("%s/dNdy_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$rundir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	$outputfile = sprintf("%s/dNdy_%s_useNucleus%s_grid%s_g2mu%s_m%s_run%05d.root",$text_out_dir,$package,$par_useNucleus,$par_size,$par_g2mu,$par_m,$irun);
	print FILE "mv outfile.root ${outputfile}\n"; 

#	print FILE "mv -v *.dat ${text_out_dir}\n";

	print FILE "rm -rf *\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";
}

