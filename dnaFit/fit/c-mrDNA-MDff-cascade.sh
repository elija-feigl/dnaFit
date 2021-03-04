############################################################################
# General input
############################################################################
# PATH
readonly UBIN=/home/gpuser/NAMD_2.13_Linux-x86_64-multicore/
readonly VBIN=/usr/local/bin/

# general
readonly DESIGNNAME=Abrickv4-2_barefob5
declare -ri TSTEPS=12000

# cascade
declare -ri NCASCADE=8
readonly GFMAX=22
readonly MAPRESOLUTION=8


############################################################################
# Advanced input
############################################################################
# location
readonly DIR=$PWD

# start with enrgMD
declare -ri TSTEPSERG=4800
declare -ri MSTEPSERG=4800
# input mrc map
readonly MAPTHRES=0.0
# mdff settings
readonly DIEL=1 #dielectric
readonly GSCALE=0.3

# repeat cascade for imporved docking
declare -ri NITER=2
declare -ri REPEATSTOP=2
declare -ri REMOVE_LR=2
# should the system equilibrate after first step
readonly RELAX=true
declare -ri TSTEPSRELAX=18000
# remove interhelical bonds after step NREFINE
readonly REFINE=true
declare -ri TSTEPSREFINE=6000
# final with increased GSCALE
readonly FINAL=true
readonly GSCALEFINAL=1.0
declare -ri MSTEPSFINAL=12000
# production run for average structure: false if final
readonly AVG=false
declare -ri TSTEPSAVG=24000

############################################################################
############################################################################
# namd run $NAMDFILE $STEP
run_namd() {
	OUTPUT=$DIR/$2
	if  [ ! -d "$OUTPUT" ]; then
		mkdir $OUTPUT
		OUT=$DIR/cMDff$2.log
		case "CUDA" in
		*$UBIN*) $UBIN/namd2 +idlepoll +p8 +devices 0,1,2,3 $1 2>&1 | tee $OUT ;;
		*)       $UBIN/charmrun +p32 $UBIN/namd2 +netpoll $1 2>&1 | tee $OUT ;;
		esac
	else
		echo "skipping step $2 (already exists)"
	fi
}

DEBUGrun_namd() {
	OUTPUT=$DIR/$2
	if  [ ! -d "$OUTPUT" ]; then
		case "CUDA" in
		*$UBIN*) echo "run CUDA" ;;
		*)       echo "run CPU" ;;
		esac
	fi
}

# change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
change_NAMDfile() {
	sed -i "s/set TS .*/set TS $1/g" $NAMDFILE
	sed -i "s/set MS .*/set MS $2/g" $NAMDFILE
	sed -i "s/set N .*/set N $3/g" $NAMDFILE
	sed -i "s/set PREVIOUS .*/set PREVIOUS $4/g" $NAMDFILE
	sed -i "s/set GRIDON .*/set GRIDON $5/g" $NAMDFILE
	sed -i "s/set GRIDFILE .*/set GRIDFILE $6/g" $NAMDFILE
	sed -i "s/set GSCALE .*/set GSCALE $7/g" $NAMDFILE
	sed -i "s/set ENRGMDON .*/set ENRGMDON $8/g" $NAMDFILE
	sed -i "s/set ENRGMDBONDS .*/set ENRGMDBONDS \$PREFIX$9/g" $NAMDFILE
	sed -i "s/set OUTPUTNAME .*/set OUTPUTNAME ${10}\/\$PREFIX/g" $NAMDFILE
}

# generate MDff grid data using VMD
generate_VMDprep() {
	# prep maps
	echo "package require volutil"  >> $VMDFILE
	echo "package require mdff"  >> $VMDFILE

	echo "volutil -clamp $MAPTHRES:1.0 $DIR/$DESIGNNAME.mrc -o $DIR/base.dx"  >> $VMDFILE
	for N in $(seq 0 $(($NCASCADE-2))); do
		if [[ "$RELAX" = true  || $N != 0 ]]; then
			GFDIFF=$((GFMAX-MAPRESOLUTION))
			GFX=$((GFDIFF/NCASCADE))
			X=$((NCASCADE-N))
			GFADD=$((GFX*X))
			GF=$((GFADD+MAPRESOLUTION))
			echo "volutil -smooth $GF $DIR/base.dx -o $DIR/$N.dx"  >> $VMDFILE
		fi
	done

	# prep mdff
	for N in $(seq 0 $(($NCASCADE-2))); do
		echo "mdff griddx -i $DIR/$N.dx -o $DIR/grid-$N.dx" >> $VMDFILE
	done
	echo "mdff griddx -i $DIR/base.dx -o $DIR/grid-$(($NCASCADE-1)).dx" >> $VMDFILE
	echo "mdff griddx -i $DIR/base.dx -o $DIR/grid-base.dx" >> $VMDFILE

	echo "mdff gridpdb -psf $DIR/$DESIGNNAME.psf -pdb $DIR/$DESIGNNAME.pdb -o $DIR/grid.pdb" >> $VMDFILE
	echo "exit" >>  $VMDFILE
}

run_vmd_prep() {
	declare VMDFILE=$DIR/mdff-prep.vmd
	touch $VMDFILE
	generate_VMDprep
	$VBIN/vmd -dispdev text -eofexit -e $VMDFILE
}

# necessary updates after each run
update_TSLAST_PREVIOUS() {
	declare -g TSLAST+=$1
	sed -i "s/set TSLAST .*/set TSLAST $TSLAST/g" $NAMDFILE
	declare -g PREVIOUS=$2
}

run_vmd_post() {
	VMDFILE2=$DIR/mdff-post.vmd
	touch $VMDFILE2

	echo "postprocessing"
	# prep packages
	echo "package require volutil"  >> $VMDFILE2
	echo "package require mdff"  >> $VMDFILE2

	# merge dcd
	echo "mol new $DESIGNNAME.psf" >> $VMDFILE2
	echo "mol addfile ./enrgMD/$DESIGNNAME.dcd start 0 step 1 waitfor all" >> $VMDFILE2
	echo "for {set i 0} {\$i < ($STEP+1)} {incr i 1} { mol addfile ./\$i/$DESIGNNAME.dcd start 0 step 1 waitfor all }" >> $VMDFILE2
	if [ "$FINAL" = true ]; then
	echo "mol addfile ./final/$DESIGNNAME.dcd start 0 step 1 waitfor all" >> $VMDFILE2
	fi
	if [ "$AVG" = true ]; then
		echo "mol addfile ./avg/$DESIGNNAME.dcd start 0 step 1 waitfor all" >> $VMDFILE2
	fi
	echo "animate write dcd $DIR/$DESIGNNAME.dcd" >> $VMDFILE2
	# last frame pdb
	echo "set sel [atomselect top all]" >> $VMDFILE2
	echo "\$sel writepdb $DIR/$DESIGNNAME-last.pdb" >> $VMDFILE2
	
	# avg pdb
	if [ "$AVG" = true ]; then
		echo "mol new $DESIGNNAME.psf" >> $VMDFILE2
		echo "mol addfile ./avg/$DESIGNNAME.dcd start 0 step 1 waitfor all" >> $VMDFILE2
		echo "set sel [atomselect top all]" >> $VMDFILE2
		echo "set pos [measure avpos \$sel]" >> $VMDFILE2
		echo "\$sel set {x y z} \$pos" >> $VMDFILE2
		echo "\$sel writepdb $DIR/$DESIGNNAME-avg.pdb" >> $VMDFILE2
	fi
	
	echo "exit" >>  $VMDFILE2

	$VBIN/vmd -dispdev text -eofexit -e $VMDFILE2
}

cleanup() {
	mkdir $DIR/data 
	mv $DIR/enrgMD $DIR/data/enrgMD
	for N in $(seq 0 $(($STEP))); do
		mv $DIR/$N $DIR/data/
	done
	if [ "$FINAL" = true ]; then
		mv $DIR/final $DIR/data/
	fi
	if [ "$AVG" = true ]; then
		mv $DIR/avg $DIR/data/
	fi

	mkdir $DIR/log
	mv $DIR/*.log $DIR/log
	mv $DIR/*.err $DIR/log
	mv $DIR/*.vmd $DIR/log

	mkdir $DIR/grid
	mv $DIR/*.dx $DIR/grid
	mv $DIR/grid.pdb $DIR/grid

	mkdir $DIR/startup
	mv $DIR/charmm36.nbfix $DIR/startup
	mv $DIR/*.namd $DIR/startup
	mv $DIR/*.sh $DIR/startup
	mv $DIR/*.exb $DIR/startup/
	mv $DIR/$DESIGNNAME.pdb $DIR/startup/
	cp $DIR/$DESIGNNAME.psf $DIR/startup/
	cp $DIR/$DESIGNNAME.mrc $DIR/startup/
}

############################################################################
############################################################################

############################################################################
#	prepare cascading maps for MDff (requires VMD with volutil, mdff)
echo "Prepare MDff"
echo $(date)

FILE=$DIR/grid.pdb
if  [ ! -f "$FILE" ]; then
	echo "prepare cascading maps for MDff"
	run_vmd_prep
else
	echo "no need to prepare cascading maps for MDff"
fi


############################################################################
#set general parameters
NAMDFILE=$DIR/CeMDff.namd
sed -i "s/set PREFIX .*/set PREFIX $DESIGNNAME/g" $NAMDFILE
sed -i "s/set DIEL .*/set DIEL $DIEL/g" $NAMDFILE
declare -i TSLAST=0
sed -i "s/set TSLAST .*/set TSLAST 0/g" $NAMDFILE

#prep intrahelical bond file # TODO: move to generation script
cp $DIR/$DESIGNNAME.exb $DIR/$DESIGNNAME-SR.exb
sed -i "/31$/d" $DIR/$DESIGNNAME-SR.exb
sed -i "/31.00$/d" $DIR/$DESIGNNAME-SR.exb

############################################################################
# energy MD step
echo "pure enrgMD relaxation"
echo $(date)
declare -i STEP=-1
CURRENT="enrgMD"
EXB=".exb"
#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
change_NAMDfile $TSTEPSERG $MSTEPSERG $STEP "0" "0" "0" "0" "on" $EXB $CURRENT
run_namd $NAMDFILE $CURRENT
update_TSLAST_PREVIOUS $TSTEPSERG $CURRENT

############################################################################
# casccade loop
echo "cascading-enrgMDff"
echo $(date)
for P in $(seq 0 $(($NITER-1))); do
	# for bad docking the system will be relxed with pure enrgMD after first iteration
	if [ "$RELAX" = true ] && [ "$NITER" -gt 1 ] && [ "$P" -eq 1 ]; then
		STEP+=1
		EXB=".exb"	
		echo "relax cascading-enrgMDff after first iteratino"
	 	#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
		change_NAMDfile $TSTEPSRELAX "0" $STEP $PREVIOUS "0" "0" "0" "on" $EXB "\$N"
		run_namd $NAMDFILE $STEP	
		update_TSLAST_PREVIOUS $TSTEPSRELAX $STEP
	fi

	if [ "$NITER" -gt 1 ] && [ "$P" -eq 0 ]; then
		NRUNS=$REPEATSTOP
	else
		NRUNS=$NCASCADE
	fi

	for N in $(seq 0 $(($NRUNS-1))); do
		STEP+=1	
		if [ "$N" -gt "$REMOVE_LR" ]; then
			EXB="-SR.exb"
		else
			EXB=".exb"
		fi
		echo "step $STEP with N $N and iteration $P"
	 	#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
		change_NAMDfile $TSTEPS "0" $STEP $PREVIOUS "1" "grid-$N.dx" $GSCALE "on" $EXB "\$N"
		run_namd $NAMDFILE $STEP	
		update_TSLAST_PREVIOUS $TSTEPS $STEP		
	done
done

############################################################################
# refine by removing intrahelical bonds
if [ "$REFINE" = true ]; then
	STEP+=1	
	echo "step $STEP with refine"
	echo $(date)
	#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
	change_NAMDfile $TSTEPSREFINE "0" $STEP $PREVIOUS "1" "grid-base.dx" $GSCALE "on" "-HO.exb" "\$N"
	run_namd $NAMDFILE $STEP
	update_TSLAST_PREVIOUS $TSTEPSREFINE $STEP
fi

############################################################################
# finalize by energy min
if [ "$FINAL" = true ]; then
	N=$((STEP+1))
	if [ "$REFINE" = true ]; then
		EXB="-HO.exb"
	else
		EXB="-SR.exb"
	fi
	echo "step $STEP with energy minimisation"
	echo $(date)
	CURRENT="final"
	#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
	change_NAMDfile "0" $MSTEPSFINAL $N $PREVIOUS "1" "grid-base.dx" $GSCALEFINAL "on" $EXB $CURRENT
	run_namd $NAMDFILE $CURRENT
	update_TSLAST_PREVIOUS $MSTEPSFINAL $CURRENT
fi

############################################################################
#average struct
if [  "$AVG" = true ] && [  "$FINAL" = false ]]; then
	N=$((STEP+2))
	CURRENT="avg"
	if [ "$REFINE" = true ]; then
		EXB="-HO.exb"
	else
		EXB="-SR.exb"
	fi
	echo "Starting production run for average structure"
	echo $(date)

	#change_NAMDfile $TS $MS $N $PREVIOUS $GRIDON $GRIDFILE $GSCALE $ENRGMDON $ENRGMDBONDS $OUTPUTNAME
	change_NAMDfile $TSTEPSAVG "0" $N $PREVIOUS "1" "grid-base.dx" $GSCALE "on" $EXB $CURRENT
	run_namd $NAMDFILE $CURRENT
	update_TSLAST_PREVIOUS $TSTEPSAVG $CURRENT
fi

############################################################################
#	merge trajectory and prepare ccc masked map
run_vmd_post

############################################################################
#	cleanup
cleanup

echo "done"
echo $(date)