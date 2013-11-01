char version[]="cmtool, v6, roberto toro, 1 July 2013";	// added labelling function based on meSH tags from brainspell
//char version[]="cmtool, v5, roberto toro, 14 June 2013";	// brainspell instead of brainmap!
//char version[]="cmtool, v4, roberto toro, 28 May 2013";	// fixed a bug in main that prevented many individual maps to be saved in abscence of the -average switch
//char version[]="cmtool, v3, roberto toro, 31 Jan 2011";	// unflip g_roi2mni (why was it flipped??), add ROIarticles, coverSeed flag, coverage threshold
//char version[]="cmtool, v2, roberto toro, 12 Dec 2010";	// accept singleton networks, fix peak detection
//"cmtool, v1, roberto toro, 21 Oct 2010";

/* TO DO
	finish seedMNI	
	how do I call a seed when the volume is averaged? (the average seed?)
	g_roi2mni matrix should be configured from defaults.txt and not hard-wired

   Notes
	Compiled using:
	cd ~/Applications/brainbits/f.cmtool
	gcc -Wall -c cmtool.c -I /Users/roberto/Documents/_03_Papers/2005_11Coactivation-CerebCortex/2007_04CoactivationMap-CerebCortex/bin/CoactivationMap -I ~/Applications/brainbits/RAMONES/
	gcc -Wall cmtool.o coactivation.o Analyze.o -o cmtool
	
	Tested commands:
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -vol 6,24,23 -convert lr -out ~/Desktop/ -save map
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -vol 6,24,23 -convert phir -out ~/Desktop/ -save map
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -roi ~/Documents/2010_06Rest-Milham/analysis/01roi/left_insula_ROI.hdr -convert phir -average -out ~/Desktop/ -save map
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -vol 6,24,23 -convert phir -out ~/Desktop/ -save peaks 0.2 4
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -vol 6,24,23 -convert phir -out ~/Desktop/ -save articles 0.2 4
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -roi ~/Documents/2010_06Rest-Milham/analysis/01roi/left_insula_ROI.hdr -convert phir -average -out ~/Desktop/ -save peaks 0.1 4
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -roi ~/Documents/2010_06Rest-Milham/analysis/01roi/left_insula_ROI.hdr -convert phir -average -coverseed -out ~/Desktop/ -save articles 0.1 4
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -roi ~/Documents/2010_06Rest-Milham/analysis/01roi/left_insula_ROI.hdr -convert phir -average -coverseed -coverage 50 -out ~/Desktop/ -save articles 0.1 4
	./cmtool /Volumes/4TB/data/coactivations/2010/coincidences -roi ~/Documents/2010_06Rest-Milham/analysis/01roi/left_insula_ROI.hdr -out ~/Desktop/ -saveroiarticles
	./cmtool ~/Documents/2010_10Coactivations/coincidences2013 -vol 6,24,23 -convert dct -out ~/Desktop/ -save map
	./cmtool ~/Documents/2010_10Coactivations/coincidences2013 -roi ~/Desktop/mask.hdr -convert dct -out ~/Desktop/cmap-dct/ -average -save map
*/

#include "defaults.h"
#include "coactivation.h"
#include "Analyze.h"

typedef struct
{
	char	tag[128];
	float	r2;
}WTag;

#define CoincidencesType			1
#define LikelihoodRatioType			2
#define	LogPValueType				3
#define	TScoreType					4
#define	MutualInformationType		5
#define	PhiCorrelationType			6
#define DiscreteCosineTransformType	7

char	*g_cmapdir;				// coactivation map directory (contains defaults.txt, sum.img, brainspell.xml/*brainmap.xml*/)
char	g_roipath[2048];		// path to ROI (used by seedROI and saveROIArticles)
char	g_out[2048];			// output files name root
short	*g_coin;				// current coincidences map
float	*g_map;					// current map (computed from g_coin)
float	*g_average;				// average map (computed by adding all g_map(g_coin) for an ROI)
short	*g_sum;					// coactivation map sum volume
int		g_nseeds;				// number of seeds
int		*g_seeds;				// 3D coordinates of seeds in volume index space	
int		g_swapFlag;				// to swap or not to swap the volumes (default 0=no)
int		g_averageFlag;			// to average or not to average the ROI volumes (default 0=no)
int		g_meshdisplay;			// display mode for mesh tags, either as "code" (1) or "name" (2)
char	g_meshroot[256];		// root used to filter the displayed mesh heading codes, default "F"
int		g_coverSeedFlag;		// in article selection, only report articles that cover the seed
float	g_coverageThreshold;	// in article selection, only report articles that cover at least this percentage of the network
char	*g_dataType[]={"coin","lr","logp","t","mi","phir","dct"};	// data type string description
int		g_dataTypeIndex;			// data type index
float	g_roi2mni[16]={ 4,0,0,0,	// conversion matrix from standard coactivation map volume (45,54,45, but should be read from defaults...) to MNI
						0,4,0,0,
						0,0,4,0,
					   -88,-124,-70,  1};
GlobalDefaults	gd;		// global defaults (to eliminate in future versions?)

void v_m(float *r,float *v,float *m);

#pragma mark -
int cmp(char *str, char *root)
{
	int	i;
	for(i=0;i<strlen(root);i++)
		if(str[i]!=root[i])
			break;
	return i;
}
void seedMNI(char *coords)
{
	float	x[3],y[3];
	
	sscanf(coords," %f , %f , %f ",&x[0],&x[1],&x[2]);
	g_nseeds=1;
	g_seeds=(int*)calloc(g_nseeds*3,sizeof(int));
	v_m(y,x,gd.t2v);
	g_seeds[0]=y[0];
	g_seeds[1]=y[1];
	g_seeds[2]=y[2];
}
#pragma mark _
void seedVol(char *coords)
{
	int	i[3];

	sscanf(coords," %i , %i , %i ",&i[0],&i[1],&i[2]);
	g_nseeds=1;
	g_seeds=(int*)calloc(g_nseeds*3,sizeof(int));
	g_seeds[0]=i[0];
	g_seeds[1]=i[1];
	g_seeds[2]=i[2];
}
#pragma mark _
float getValue(AnalyzeHeader *hdr, char *img, int x, int y, int z)
{
	float		val;
	RGBValue	rgb;

	if(hdr->datatype==RGB)
	{
		rgb=((RGBValue*)img)[z*hdr->dim[2]*hdr->dim[1]+y*hdr->dim[1]+x];
		val=((int)rgb.r)>>16|((int)rgb.g)>>8|((int)rgb.b);
	}
	else
	{
		switch(hdr->datatype)
		{	case UCHAR: val=((unsigned char*)img)[z*hdr->dim[2]*hdr->dim[1]+y*hdr->dim[1]+x];	break;
			case SHORT: val=        ((short*)img)[z*hdr->dim[2]*hdr->dim[1]+y*hdr->dim[1]+x];	break;
			case INT:	val=          ((int*)img)[z*hdr->dim[2]*hdr->dim[1]+y*hdr->dim[1]+x];	break;
			case FLOAT:	val=        ((float*)img)[z*hdr->dim[2]*hdr->dim[1]+y*hdr->dim[1]+x];	break;
		}
	}
	return val;
}
void seedROI(char *path)
{
	printf(" > seedROI\n");
	char			*addr;
	AnalyzeHeader	*hdr;
	char			*img;
	int				sz,swapped;
	int				i,j,k;
	float			mni[3],voxroi[3],voxcmap[3];
	
	Analyze_load(path,&addr,&sz,&swapped);
	hdr=(AnalyzeHeader*)addr;
	img=(char*)(addr+sizeof(AnalyzeHeader));
	
	g_nseeds=0;
	for(i=0;i<hdr->dim[1];i++)
	for(j=0;j<hdr->dim[2];j++)
	for(k=0;k<hdr->dim[3];k++)
	if(getValue(hdr,img,i,j,k)>0)
		g_nseeds++;
	g_seeds=(int*)calloc(g_nseeds*3,sizeof(int));
	printf("%i seeds\n",g_nseeds);

	/* // Print transformation matrices
	for(j=0;j<4;j++){for(i=0;i<4;i++) printf("%g\t",g_roi2mni[j*4+i]); printf("\n");}
	printf("\n");
	for(j=0;j<4;j++){for(i=0;i<4;i++) printf("%g\t",gd.t2v[j*4+i]);printf("\n");}
	*/

	g_nseeds=0;
	for(i=0;i<hdr->dim[1];i++)
	for(j=0;j<hdr->dim[2];j++)
	for(k=0;k<hdr->dim[3];k++)
	if(getValue(hdr,img,i,j,k)>0)
	{
		voxroi[0]=i;
		voxroi[1]=j;
		voxroi[2]=k;
		v_m(mni,voxroi,g_roi2mni);
		v_m(voxcmap,mni,gd.t2v);
		g_seeds[g_nseeds*3+0]=voxcmap[0];
		g_seeds[g_nseeds*3+1]=voxcmap[1];
		g_seeds[g_nseeds*3+2]=voxcmap[2];
		g_nseeds++;
		
		// Print conversion: // printf("%g,%g,%g -> %g,%g,%g -> %g,%g,%g\n",voxroi[0],voxroi[1],voxroi[2],mni[0],mni[1],mni[2],voxcmap[0],voxcmap[1],voxcmap[2]);
	}
	free(addr);
}
#pragma mark -
void addToAverage(void)
{
	int		sz,i;

	sz=gd.LR*gd.PA*gd.IS;
	if(g_average==NULL)
		g_average=(float*)calloc(sz,sizeof(float));
	for(i=0;i<sz;i++)
		g_average[i]+=g_map[i];
}
void average(void)
{
	int		sz,i;

	sz=gd.LR*gd.PA*gd.IS;
	for(i=0;i<sz;i++)
		g_map[i]=g_average[i]/(float)g_nseeds;
	free(g_average);
	g_average=NULL;
}

#pragma mark -
void convert(int j)
{
	printf(" > convert [%s]\n",g_dataType[g_dataTypeIndex-1]);
	int		i1,i2,sz;

	i1=g_seeds[3*j+2]*gd.PA*gd.LR+g_seeds[3*j+1]*gd.LR+g_seeds[3*j+0];
	sz=gd.LR*gd.PA*gd.IS;
	for(i2=0;i2<sz;i2++)
	{
		switch(g_dataTypeIndex)
		{
			case LikelihoodRatioType:
				g_map[i2]=likelihood_ratio(g_sum[i1], g_sum[i2], g_coin[i2], gd.N);
				break;
			case LogPValueType:
				g_map[i2]=-log10(p_combination(g_sum[i1], g_sum[i2], g_coin[i2], gd.N));
				break;
			case TScoreType:
				g_map[i2]=t_score(g_sum[i1], g_sum[i2], g_coin[i2], gd.N);
				break;
			case MutualInformationType:
				g_map[i2]=mutual_information(g_sum[i1], g_sum[i2], g_coin[i2], gd.N);
				break;
			case PhiCorrelationType:
			case DiscreteCosineTransformType:
				g_map[i2]=phi_correlation(g_sum[i1], g_sum[i2], g_coin[i2], gd.N);
				break;
		}
	}
	
	if(g_dataTypeIndex==DiscreteCosineTransformType)
	{
		int	dim[3];
		dim[0]=gd.LR;
		dim[1]=gd.PA;
		dim[2]=gd.IS;
		discrete_cosine_transform(g_map,dim);
	}
}
#pragma mark -
void saveMap(int i, int j, int k)
{
	printf(" > saveMap\n");
	char	path[2048];
	FILE	*f;
	int		sz;
	
	sz=gd.LR*gd.PA*gd.IS;
	if(g_averageFlag)
		sprintf(path,"%saverage-%s.img",g_out,g_dataType[g_dataTypeIndex-1]);
	else
		sprintf(path,"%s%03i%03i%03i-%s.img",g_out,i,j,k,g_dataType[g_dataTypeIndex-1]);

	printf("%s\n",path);
	f=fopen(path,"w");
	fwrite(g_map,sz,sizeof(float),f);
	fclose(f);
}
#pragma mark _
void savePeaks(int npeaks, Peak *peaks, int j)
{
	printf(" > savePeaks\n");
	int		i,i2;
	FILE	*f;
	char	path[2048];
	float	x[3],y[3];
	
	sprintf(path,"%speaks%i.txt",g_out,j);
	f=fopen(path,"w");
	fprintf(f,"I\tJ\tK\tX\tY\tZ\tMaxExp\tMaxStat\n");
	for(i=0;i<npeaks;i++)
	{
		x[0]=peaks[i].a;
		x[1]=peaks[i].b;
		x[2]=peaks[i].c;
		i2=x[2]*gd.PA*gd.LR+x[1]*gd.LR+x[0];
		v_m(y,x,gd.v2t);
		fprintf(f,"%i\t%i\t%i\t%g\t%g\t%g\t%i\t%g\n",
					peaks[i].a,peaks[i].b,peaks[i].c,
					y[0],y[1],y[2],
					g_coin[i2],g_map[i2]);
	}
	fclose(f);
}
#pragma mark _
void configureArticlesLUT(int *LUT)
{
	FILE	*f;
	char	path[2048],str[2048];
	int		papPos,nexp;
	int		flagExperiment;
	
	sprintf(path,"%s/brainspell.xml",g_cmapdir);
	//sprintf(path,"%s/brainmap.xml",g_cmapdir);
	f=fopen(path,"r");
	flagExperiment=0;
	nexp=0;
	while(!feof(f))
	{
		fgets(str,2048,f);
		
		if(strstr(str,"<paper>"))
			papPos=ftell(f);
		else
		if(strstr(str,"<experiment>"))
		{
			LUT[2*nexp+0]=papPos;		// beginning of the paper
			LUT[2*nexp+1]=ftell(f);		// beginning of the experiment
			nexp++;
		}
	}
	fclose(f);
}
void articleDescription(int *LUT,int ind,char *pubmedid,char *reference,char *meshcodes,char *domains)
{
	FILE	*f;
	char	path[2048],str[2048],*ptr0,*ptr1;
	int		nexp,firstAuthorFlag;
	char	author[256],year[256],domain[256],code[256],name[256];
	int		flagExperiment,flagFinished;
	
	sprintf(path,"%s/brainspell.xml",g_cmapdir);
	//sprintf(path,"%s/brainmap.xml",g_cmapdir);
	f=fopen(path,"r");
	
	// get paper info
	fseek(f,LUT[2*ind+0],SEEK_SET);
	firstAuthorFlag=0;
	flagFinished=0;
	strcpy(meshcodes,"");
	do
	{
		fgets(str,2048,f);
		if(strstr(str,"<Medline_number>"))
		{
			ptr0=strstr(str,"<Medline_number>");
			strcpy(pubmedid,ptr0+strlen("<Medline_number>"));
			ptr1=strstr(pubmedid,"</Medline_number>");
			ptr1[0]=(char)0;
		}
		else
		if(firstAuthorFlag==0 && strstr(str,"<author>"))
		{
			ptr0=strstr(str,"<author>");
			strcpy(author,ptr0+strlen("<author>"));
			ptr1=strstr(author,"</author>");
			ptr1[0]=(char)0;
			firstAuthorFlag=1;
		}
		else
		if(strstr(str,"<year>"))
		{
			ptr0=strstr(str,"<year>");
			strcpy(year,ptr0+strlen("<year>"));
			ptr1=strstr(year,"</year>");
			ptr1[0]=(char)0;
		}
		else
		if(strstr(str,"<name>"))
		{
			ptr0=strstr(str,"<name>");
			strcpy(name,ptr0+strlen("<name>"));
			ptr1=strstr(name,"</name>");
			ptr1[0]=(char)0;
		}
		else
		if(strstr(str,"<code>"))
		{
			ptr0=strstr(str,"<code>");
			strcpy(code,ptr0+strlen("<code>"));
			ptr1=strstr(code,"</code>");
			ptr1[0]=(char)0;
			if(cmp(code,g_meshroot)==strlen(g_meshroot))
			{
				if(meshcodes[0]==(char)0)
					strcpy(meshcodes,(g_meshdisplay==1)?code:name);
				else
					sprintf(meshcodes,"%s;%s",meshcodes,(g_meshdisplay==1)?code:name);
			}
		}
		else
		if(strstr(str,"</MeshCodes>"))
		{
			flagFinished=1;
		}
	}
	while(flagFinished==0 && !feof(f));
	sprintf(reference,"\"%s, %s\"",author,year);

	// get behavioural domain info
	strcpy(domains,"");
	fseek(f,LUT[2*ind+1],SEEK_SET);
	flagExperiment=1;
	nexp=0;
	while(!feof(f) && flagExperiment==1)
	{
		fgets(str,2048,f);
		
		if(strstr(str,"<domain>"))
		{
			ptr0=strstr(str,"<domain>");
			strcpy(domain,ptr0+strlen("<domain>"));
			ptr1=strstr(domain,"</domain>");
			ptr1[0]=(char)0;
			if(domains[0]==(char)0)
				strcpy(domains,domain);
			else
				sprintf(domains,"%s,%s",domains,domain);
		}
		else
		if(strstr(str,"</experiment>"))
			flagExperiment=0;
	}
	fclose(f);
}
void saveArticles(int npeaks, Peak *peaks, int j)
{
	printf(" > saveArticles (coverseed=%s)\n",g_coverSeedFlag?"YES":"NO");
	char	*xvol;
	int		i,k,l,sz,npks;
	char	path[2048];
	char	name[255];
	FILE	*f,*froi;
	int		*LUT;
	char	pubmedid[256],reference[256],domains[2048],meshcodes[2048];
	
	LUT=(int*)calloc(2*gd.N, sizeof(int));
	configureArticlesLUT(LUT);
	
	sz=gd.LR*gd.PA*gd.IS;
	
	sprintf(path,"%sarticles%i.txt",g_out,j);
	f=fopen(path,"w");
	fprintf(f,"Peaks\tPercent\tExpNumber\tPubMedID\tReference\tMeshCodes\tDomains\n");
	// scan through the experiments for intersections with the network peaks
	xvol=(char*)calloc(sz,sizeof(char));
	for(i=1;i<=gd.N;i++)
	{
		if((i%(gd.N/100))<((i-1)%(gd.N/100)))
		{
			printf("%i%% ",(int)(100*i/(float)gd.N));
			fflush(stdout);
		}
		sprintf(name,"%s/rois/%i.img",g_cmapdir,i);
		froi=fopen(name,"r");
		fread(xvol,sz,sizeof(char),froi);
		fclose(froi);
		
		npks=0;
		for(k=0;k<npeaks;k++)
		if(xvol[peaks[k].c*gd.PA*gd.LR+peaks[k].b*gd.LR+peaks[k].a])
			npks++;
		
		// check seed covering
		if(g_coverSeedFlag)
		{
			int	cover=0;
			for(l=0;l<g_nseeds;l++)
				cover+=xvol[g_seeds[3*l+2]*gd.PA*gd.LR+g_seeds[3*l+1]*gd.LR+g_seeds[3*l+0]];
			if(cover==0)
				continue;
		}
		
		// check minimum covering
		if(npks>=1 && 100*npks/(float)npeaks>g_coverageThreshold)
		{
			articleDescription(LUT,i-1,pubmedid,reference,meshcodes,domains);
			fprintf(f,"%i\t%g\t%i\t%s\t%s\t%s\t%s\n",npks,100*npks/(float)npeaks,i,pubmedid,reference,meshcodes,domains);
		}
	}
	fclose(f);
	free(LUT);
}
int compare(const void *a, const void *b)
{
	float	va,vb;
	va=((WTag*)a)->r2;
	vb=((WTag*)b)->r2;
	
	if(va>vb)
		return -1;
	if(va<vb)
		return 1;
	return 0;
}
void saveTopMesh(int jj)
{
	printf(" > saveTopMesh\n");
	short	*xvol;
	int		i,l,sz;
	char	path[2048],cmd[2048];
	FILE	*f,*ftag,*fmesh;
	int		ntag;
	char	tag[512];
	float	x,xx,y,yy,xy,n,num,den,r2;
	WTag	*wtag;
	char	*tmpdir;
	
	sz=gd.LR*gd.PA*gd.IS;

	// total map mass
	x=xx=0;
	for(i=0;i<gd.LR*gd.PA*gd.IS;i++)
	{
		x+=g_map[i];
		xx+=g_map[i]*g_map[i];
	}
	n=gd.LR*gd.PA*gd.IS;

	// scan through the experiments
	wtag=(WTag*)calloc(1,sizeof(WTag));
	xvol=(short*)calloc(sz,sizeof(short));
	sprintf(path,"%s/mesh.txt",g_cmapdir);
	fmesh=fopen(path,"r");

	ntag=0;
	while(!feof(fmesh))
	{
		fgets(tag,512,fmesh);
		ntag++;
	}
	fclose(fmesh);

	wtag=(WTag*)calloc(ntag,sizeof(WTag));	
	fmesh=fopen(path,"r");
	tmpdir=tmpnam(NULL);
	sprintf(cmd,"mkdir %s",tmpdir);
	system(cmd);
	for(l=0;l<ntag;l++)
	{
		if((l%(ntag/100))<((l-1)%(ntag/100)))
		{
			printf("%i%% ",(int)(100*l/(float)ntag));
			fflush(stdout);
		}

		fgets(tag,512,fmesh);
		if(tag[strlen(tag)-1]=='\r'||tag[strlen(tag)-1]=='\n')
			tag[strlen(tag)-1]=(char)0;
		
		// open tag file
		sprintf(path,"%s/tags.zip",g_cmapdir);
		sprintf(cmd,"unzip -joq %s \"%s.img\" -d %s/",path,tag,tmpdir);
		system(cmd);
		sprintf(path,"%s/%s.img",tmpdir,tag);
		ftag=fopen(path,"r");
		fread(xvol,sz,sizeof(short),ftag);
		fclose(ftag);
		
		// compute correlation
		y=yy=xy=0;
		for(i=0;i<sz;i++)
		{
			y+=xvol[i];
			yy+=xvol[i]*xvol[i];
			xy+=g_map[i]*xvol[i];
		}
		num=xy-x*y/n;
		den=sqrt((xx-x*x/n)*(yy-y*y/n));
		if(den==0)
			r2=0;
		else
			r2=pow(num/den,2);
		
		strcpy(wtag[l].tag,tag);
		wtag[l].r2=r2;
	}
	printf("\n");
	sprintf(cmd,"rm -r %s",tmpdir);
	system(cmd);
	
	// print top 100
	qsort(wtag,ntag,sizeof(WTag),compare);
	sprintf(path,"%stop100.%i.txt",g_out,jj);
	f=fopen(path,"w");
	fprintf(f,"R2\tMeshCode\n");
	for(l=0;l<100;l++)
		fprintf(f,"%f\t%s\n",wtag[l].r2,wtag[l].tag);
	fclose(f);
	free(wtag);
	
}
#pragma mark -
void saveROIArticles(void)
{
	printf(" > saveROIArticles\n");
	char			path[2048];
	char			*xvol;
	int				i,n,sz,nvox;
	char			name[255];
	FILE			*f,*froi;
	int				*LUT;
	char			pubmedid[256],reference[256],domains[2048],meshcodes[2048];
	
	LUT=(int*)calloc(2*gd.N, sizeof(int));
	configureArticlesLUT(LUT);
	
	sz=gd.LR*gd.PA*gd.IS;
	
	sprintf(path,"%sroiarticles.txt",g_out);
	f=fopen(path,"w");
	fprintf(f,"Voxels\tExpNumber\tPubMedID\tReference\tMeshCodes\tDomains\n");
	// scan through the experiments for intersections with the ROI non-null voxels
	xvol=(char*)calloc(sz,sizeof(char));
	for(n=1;n<=gd.N;n++)
	{
		if((n%(gd.N/100))<((n-1)%(gd.N/100)))
		{
			printf("%i%% ",(int)(100*n/(float)gd.N));
			fflush(stdout);
		}
		sprintf(name,"%s/rois/%i.img",g_cmapdir,n);
		froi=fopen(name,"r");
		fread(xvol,sz,sizeof(char),froi);
		fclose(froi);
		
		nvox=0;
		for(i=0;i<g_nseeds;i++)
			nvox+=xvol[g_seeds[3*i+2]*gd.PA*gd.LR+g_seeds[3*i+1]*gd.LR+g_seeds[3*i+0]];
		
		if(nvox>=1)
		{
			articleDescription(LUT,n-1,pubmedid,reference,meshcodes,domains);
			fprintf(f,"%i\t%i\t%s\t%s\t%s\t%s\n",nvox,n,pubmedid,reference,meshcodes,domains);
		}
	}
	printf("\n");
	free(LUT);
}
#pragma mark -
void saveTagVolume(char *theTag)
{
	printf(" > getTag: [%s]\n",theTag);
	char	*xvol;
	short	*sumVol;
	int		i,l,sz;
	char	path[2048];
	char	name[255];
	FILE	*f,*froi;
	int		*LUT;
	char	pubmedid[256],reference[256],domains[2048],meshcodes[2048],*tag;
	
	LUT=(int*)calloc(2*gd.N, sizeof(int));
	configureArticlesLUT(LUT);
	
	sz=gd.LR*gd.PA*gd.IS;
	
	xvol=(char*)calloc(sz,sizeof(char));
	sumVol=(short*)calloc(sz,sizeof(short));
	for(l=1;l<=gd.N;l++)
	{
		if(0)
		if((l%(gd.N/100))<((l-1)%(gd.N/100)))
		{
			printf("%i%% ",(int)(100*l/(float)gd.N));
			fflush(stdout);
		}
		sprintf(name,"%s/rois/%i.img",g_cmapdir,l);
		froi=fopen(name,"r");
		fread(xvol,sz,sizeof(char),froi);
		fclose(froi);
		
		articleDescription(LUT,l-1,pubmedid,reference,meshcodes,domains);
		tag=strtok(meshcodes,";");
		while(tag!=NULL)
		{
			if(strcmp(tag,theTag)==0)
			{
				printf("%i,%s ",l,pubmedid);
				for(i=0;i<sz;i++)
					sumVol[i]+=xvol[i];
			}
			tag=strtok(NULL,";");
		}
	}
	printf("\n");
	free(LUT);
	free(xvol);
	
	float	x,xx,y,yy,xy,n=sz;
	x=xx=y=yy=xy=0;
	for(i=0;i<sz;i++)
	{
		x+=g_map[i];
		xx+=g_map[i]*g_map[i];
		y+=sumVol[i];
		yy+=sumVol[i]*sumVol[i];
		xy+=g_map[i]*sumVol[i];
	}
	printf("x:%f, xx:%f, y:%f, yy:%f, xy:%f\n",x,xx,y,yy,xy);
	printf("r=%f\n",(xy-x*y/n)/sqrt((xx-x*x/n)*(yy-y*y/n)));

	sprintf(path,"%s%s.img",g_out,theTag);
	f=fopen(path,"w");
	fwrite(sumVol,sz,sizeof(short),f);
	fclose(f);
	free(sumVol);
}
#pragma mark -
// result_vector = vector x matrix
void v_m(float *r,float *v,float *m)
{
	// v=1x3
	// m=4x4
	// r=1x3
	r[0]=v[0]*m[0*4+0]+v[1]*m[1*4+0]+v[2]*m[2*4+0] + m[3*4+0];
	r[1]=v[0]*m[0*4+1]+v[1]*m[1*4+1]+v[2]*m[2*4+1] + m[3*4+1];
	r[2]=v[0]*m[0*4+2]+v[1]*m[1*4+2]+v[2]*m[2*4+2] + m[3*4+2];
}
void findPeaks(float threshold,int R, int *npeaks, Peak *peaks)
{
	printf(" > findPeaks: ");
	int		sz;
	int		i,j,k,l,m,n,nn,i2;
	Peak	coord;
	float	val,max;
	
	sz=gd.LR*gd.PA*gd.IS;
	*npeaks=0;

	for(i=0;i<gd.LR;i++)
	for(j=0;j<gd.PA;j++)
	for(k=0;k<gd.IS;k++)
	{
		i2=k*gd.PA*gd.LR+j*gd.LR+i;
		val=g_map[i2];
		if(val>threshold)
		{
			nn=0;
			
			max=0;
			for(l=-R;l<=R;l++)
			for(m=-R;m<=R;m++)
			for(n=-R;n<=R;n++)
			if(l*l+m*m+n*n<=R*R)
				if(i+l>=0 && i+l<gd.LR &&
				   j+m>=0 && j+m<gd.PA &&
				   k+n>=0 && k+n<gd.IS)
				{
					i2=(k+n)*gd.PA*gd.LR+(j+m)*gd.LR+(i+l);
					val=g_map[i2];
					
					if(	val>threshold &&	// count immediate superthreshold neighbours
						fabs(l)<2 && fabs(m)<2 && fabs(n)<2)
						nn++;
					
					if(val>=max)
					{
						max=val;
						coord=(Peak){i+l,j+m,k+n,0,0};
					}
				}
			
			i2=k*gd.PA*gd.LR+j*gd.LR+i;
			val=g_map[i2];
			if( nn>1 && // filter single voxel clusters
				val==max &&
				coord.a==i && coord.b==j && coord.c==k)
			{
				peaks[*npeaks].a=i;
				peaks[*npeaks].b=j;
				peaks[*npeaks].c=k;
				peaks[*npeaks].maxlr=max;
				peaks[*npeaks].maxk=g_coin[i2];
				(*npeaks)++;
			}
		}
	}
	printf("%i found\n",*npeaks);
}
void swapShort(short *vol, int n)
{
	int	i;
	char	w[2];

	for(i=0;i<n;i++)
	{
		w[0]=((char*)&vol[i])[1];
		w[1]=((char*)&vol[i])[0];
		vol[i]=*(short*)w;
	}
}
void loadSum(void)
{
	FILE	*f;
	char	path[2048];
	
	sprintf(path,"%s/sum.img",g_cmapdir);
	
	// Load sum
	f=fopen(path,"r");
	fread(g_sum,gd.LR*gd.PA*gd.IS,sizeof(short),f);
	fclose(f);
	if(g_swapFlag)
		swapShort(g_sum,gd.LR*gd.PA*gd.IS);
}
void loadCoin(int j)
{
	char	cmd[2048],path[2048];
	FILE	*f;
	int		sz;
	
	sz=gd.LR*gd.PA*gd.IS;
	sprintf(path,"%s/coincidences.zip",g_cmapdir);
	sprintf(cmd,"unzip -joq %s %03i%03i%03i.img -d tmp/",path,g_seeds[3*j+0],g_seeds[3*j+1],g_seeds[3*j+2]);
	system(cmd);
	sprintf(path,"tmp/%03i%03i%03i.img",g_seeds[3*j+0],g_seeds[3*j+1],g_seeds[3*j+2]);
	f=fopen(path,"r");
	if(f==0)
	{
		printf("ERROR: file [%s] not found\n",path);
		exit(1);
	}
	fread(g_coin,sz,sizeof(short),f);
	fclose(f);
	if(g_swapFlag)
		swapShort(g_coin,gd.LR*gd.PA*gd.IS);
}
#pragma mark -
int main(int argc, char *argv[])
{
	printf("%s\n",version);

	// argv[1]: cmapdir										| cmapdir=coactivation map directory (contains defaults.txt, sum.img, brainspell.xml/*brainmap.xml*/)
	// argv[2,3]:{mni x,y,z|vol i,j,k|roi path|tag name}	| mni x,y,z= seed coordinates in MNI space,
	//														| vol i,j,k= seed coordinates in volume index space,
	//														| roi {average} path= path to ROI file to use as multiple seeds or average through all seeds
	//														| tag name= name of the tag to pull
	// argv[4...] commands:
	//		-swap											| swap endianness in coactivation map files (coincidence, sum)
	//		-average										| average over the ROI
	//		-convert {lr|mi|logp|t|phir|dct}				| convert map to lr=likelihood ratio, mi=mutual information,
	//														|    logp=log p-value of likelihood ratio test, t=t-value, phir=phi-correlation,
	//														|	 dct=discrete cosinus transform
	//		-out root										| root=root name for output files (end with / for directory)
	//		-meshdisplay {name|code}						| display mesh tags as names or codes (default is "codes")
	//		-meshfilter	root								| in article selection, only report mesh heading codes (not names) under root (default "F")
	//		-coverseed										| in article selection, only report articles that cover the seed node (default, NO)
	//		-coverage thrs									| in article selection, minimum coverage (from 0 to 100) of network nodes to report an article (default, 0)
	//		-save map										| save the map
	//	[DEL]	-save average								| save the average map (only for argv[2]=roi)
	//		-save peaks thrs R								| save list of map peaks detected with threshold thrs and radius R
	//		-save articles thrs R							| save list of articles intersecting the map peaks detected with threshold thrs and radius R
	//		-save topmesh									| rank mesh tags based on a map-weighted sum throughout all papers. Save the top 100
	//		-saveroiarticles								| save list of articles intersecting the ROI (arg[2,3]=roi path)

	int		i,j;
	char	path[2048];

	g_nseeds=0;
	g_swapFlag=0;
	g_averageFlag=0;
	g_meshdisplay=1;
	strcpy(g_meshroot,"F");
	g_coverSeedFlag=0;
	g_coverageThreshold=0;
	g_sum=NULL;
	g_average=NULL;

	// 1. configure defaults
	g_cmapdir=argv[1];
	sprintf(path,"%s/defaults.txt",g_cmapdir);
	defaults(path,&gd);
	strcpy(g_out,"./");
	
	// 3. allocate memory
	g_sum=(short*)calloc(gd.LR*gd.PA*gd.IS,sizeof(short));
	g_coin=(short*)calloc(gd.LR*gd.PA*gd.IS,sizeof(short));
	g_map=(float*)calloc(gd.LR*gd.PA*gd.IS,sizeof(float));
	
	loadSum();
		
	for(i=2;i<argc;i++)
	{
		printf("[ %i. %s\n",i,argv[i]);
		
		if(strcmp(argv[i],"-mni")==0)
			seedMNI(argv[++i]);
		else
		if(strcmp(argv[i],"-vol")==0)
			seedVol(argv[++i]);
		else
		if(strcmp(argv[i],"-roi")==0)
		{
			strcpy(g_roipath,argv[++i]);
			seedROI(g_roipath);
		}
		else
		if(strcmp(argv[i],"-tag")==0)
			saveTagVolume(argv[++i]);
		else
		if(strcmp(argv[i],"-swap")==0)
			g_swapFlag=1;
		else
		if(strcmp(argv[i],"-average")==0)
			g_averageFlag=1;
		else
		if(strcmp(argv[i],"-convert")==0)
		{
			if(strcmp(argv[i+1],"lr")==0)
				g_dataTypeIndex=LikelihoodRatioType;
			else
			if(strcmp(argv[i+1],"logp")==0)
				g_dataTypeIndex=LogPValueType;
			else
			if(strcmp(argv[i+1],"t")==0)
				g_dataTypeIndex=TScoreType;
			else
			if(strcmp(argv[i+1],"mi")==0)
				g_dataTypeIndex=MutualInformationType;
			else
			if(strcmp(argv[i+1],"phir")==0)
				g_dataTypeIndex=PhiCorrelationType;
			else
			if(strcmp(argv[i+1],"dct")==0)
				g_dataTypeIndex=DiscreteCosineTransformType;
			else
			{
				printf("ERROR: Unknown value type \"%s\"\n",argv[3]);
				return 1;
			}
			i+=1;
		}
		else
		if(strcmp(argv[i],"-out")==0)
		{
			strcpy(g_out,argv[i+1]);
			i+=1;
		}
		else
		if(strcmp(argv[i],"-meshdisplay")==0)
		{
			if(strcmp(argv[i+1],"code")==0)
				g_meshdisplay=1;
			else
			if(strcmp(argv[i+1],"name")==0)
				g_meshdisplay=2;
			else
			{
				printf("ERROR: Unknown meshdisplay type \"%s\"\n",argv[i+1]);
				return 1;
			}
			i+=1;
		}
		else
		if(strcmp(argv[i],"-meshfilter")==0)
		{
			strcpy(g_meshroot,argv[i+1]);
			i+=1;
		}
		else
		if(strcmp(argv[i],"-coverseed")==0)
			g_coverSeedFlag=1;
		else
		if(strcmp(argv[i],"-coverage")==0)
		{
			g_coverageThreshold=atof(argv[i+1]);
			i+=1;
		}
		else
		if(strcmp(argv[i],"-saveroiarticles")==0)
			saveROIArticles();
		else
		if(strcmp(argv[i],"-save")==0)
		{
			int saveType=0;
			if(strcmp(argv[i+1],"map")==0)
				saveType=1;
			else
			if(strcmp(argv[i+1],"peaks")==0)
				saveType=2;
			else
			if(strcmp(argv[i+1],"articles")==0)
				saveType=3;
			else
			if(strcmp(argv[i+1],"topmesh")==0)
				saveType=4;

			for(j=0;j<g_nseeds;j++)
			{
				loadCoin(j);
				convert(j);

				if(g_averageFlag)
				{
					addToAverage();
					if(j==g_nseeds-1)
						average();
				}
				
				if(g_averageFlag==0 || j==g_nseeds-1)
				{
					switch(saveType)
					{
						case 1: //map
						{
							saveMap(g_seeds[3*j+0],g_seeds[3*j+1],g_seeds[3*j+2]);
							break;
						}
						case 2: //peaks
						{
							float	thrs;
							int		R;
							int		npeaks;
							Peak	peaks[1000];
		
							thrs=atof(argv[i+2]);
							R=atoi(argv[i+3]);
							findPeaks(thrs,R,&npeaks,peaks);
							savePeaks(npeaks,peaks,j);
							break;
						}
						case 3: //articles
						{
							float	thrs, R;
							int		npeaks;
							Peak	peaks[1000];
		
							thrs=atof(argv[i+2]);
							R=atof(argv[i+3]);
							findPeaks(thrs,R,&npeaks,peaks);
							saveArticles(npeaks,peaks,j);
							break;
						}
						case 4: //topmesh
						{
							saveTopMesh(j);
							break;
						}
						default:
						{
							printf("ERROR: Unrecognised data type for saving (neither articles, map nor peaks)\n");
							return 1;
						}
					}
				}
			}
			
			switch(saveType)
			{
				case 1: i+=1; break;
				case 2: i+=3; break;
				case 3: i+=3; break;
				case 4: i+=1; break;
			}
		}
		else
		{
			printf("ERROR: Unknown switch \"%s\"\n",argv[i]);
			return 1;
		}
	}
	
	free(g_sum);
	free(g_coin);
	free(g_map);
		
	return 0;
}
