#!/usr/bin/perl -w

use Cwd;

$maindir = getcwd();
$dataset = "pPb_8160GeV";
$grp = 5;

$grpdir = sprintf("%s/%s_3_grp%03d",$maindir,$dataset,$grp);
mkdir $grpdir;
chdir $grpdir;

for ($iseg=0; $iseg<100; $iseg++){

#	sleep 5;

	$wrkdir = sprintf("%s/wrk_%04d",$grpdir,$iseg);
	mkdir $wrkdir;
	chdir $wrkdir;

#	symlink "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/condor_phenix", "condor";
	symlink "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/condor", "condor";
#	symlink "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/condor_highmem", "condor";

	open(FILE,">jobscript");
	print FILE "#!/usr/bin/csh -f\n\n";

#	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/sphenix/setup.csh\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/setup.csh\n";
	print FILE "echo \"PATH=\$PATH\"\n";
	print FILE "echo \"LD_LIBRARY_PATH=\$LD_LIBRARY_PATH\"\n\n";

	$condor_dir = sprintf("\$_CONDOR_SCRATCH_DIR/%s",$wrkdir);
	print FILE "mkdir -p $condor_dir\n";
	print FILE "cd $condor_dir\n";
#	print FILE "cd $wrkdir\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/source/sonicO2/template/* .\n";
	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/source/sonic/template/* .\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/source/sonic_debug/template/* .\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_dAu_grp00${grp}/event${iseg}.dat input/inited.dat\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_dAu_grp00${grp}/event${iseg}.root quarkdist.root\n\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_dAu_grp0${grp}/event${iseg}.dat input/inited.dat\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_dAu_grp0${grp}/event${iseg}.root quarkdist.root\n\n";

	print FILE "echo GENERATE\n";
	print FILE "./generate | tee generate.log\n\n";

	print FILE "echo INITE\n";
	print FILE "./initE | tee initE.log\n\n";

#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_%s_2_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_%s_2_%02d/event%d.root",$dataset,$grp,$iseg);
#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_1_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_1_%02d/event%d.root",$dataset,$grp,$iseg);
#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_%s_2_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/initedFiles_%s_2_%02d/event%d.root",$dataset,$grp,$iseg);
#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_100_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_100_%02d/event%d.root",$dataset,$grp,$iseg);
#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_0_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_0_%02d/event%d.root",$dataset,$grp,$iseg);
#	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_0_%02d/event%d.dat",$dataset,$grp,$iseg);
#	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_0_%02d/event%d.root",$dataset,$grp,$iseg);
	$filename0 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_1_%02d/event%d.dat",$dataset,$grp,$iseg);
	$filename1 = sprintf("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/code/gen_4He_O/initedFiles_%s_1_%02d/event%d.root",$dataset,$grp,$iseg);

	print FILE "cp -av $filename0 input/inited.dat\n";
	print FILE "cp -av $filename1 quarkdist.root\n\n";

	print FILE "echo VH2\n";
	print FILE "./vh2 | tee vh2.log\n\n";

#	print FILE "echo B3D\n";
#	print FILE "./b3d default > b3d.log\n\n";

#	print FILE "echo ANALYZE\n";
#	print FILE "./analyze default > ana.log\n\n";

	print FILE "cp -av analysis/default/cent0to5/details/* $wrkdir\n";
	print FILE "cp -av *.log $wrkdir\n";
	print FILE "cp -av data/snapshot $wrkdir\n";
	print FILE "rm -rf *\n\n";

	close(FILE);
	chmod 0755, "jobscript";

	system "condor_submit condor";

}
