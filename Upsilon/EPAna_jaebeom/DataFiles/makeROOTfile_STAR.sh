g++ -o runConvSTARNpart `root-config --cflags --evelibs` convert_to_root_STAR_npart.cc
g++ -o runConvSTARPt `root-config --cflags --evelibs` convert_to_root_STAR_pT.cc
g++ -o runConvCMS276Npart `root-config --cflags --evelibs` convert_to_root_CMS_npart_2p76TeV.cc
./runConvSTARNpart
./runConvSTARPt
./runConvCMS276Npart
