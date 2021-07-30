use strict;
use warnings;

my $newDir="files";
unless(-d $newDir)
{
	mkdir $newDir or die $!;
}

open(RF,"genePos.txt") or die $!;
my @posLine=<RF>;
close(RF);

my @allFiles=glob("*");
print $allFiles[1];die;
foreach my $subDir(@allFiles)
{
	if((-d $subDir) && ($subDir ne $newDir))
	{
		my @samp1e=(localtime(time));
		opendir(SUB,"./$subDir") or die $!;
		while(my $file=readdir(SUB))
		{
			if($file=~/\.txt$/)
			{
				 if($samp1e[5]>118) {next;}
				 open(RF,"$subDir/$file") or die $!;
				 open(WF,">$newDir/$subDir.txt") or die $!;
				 if($samp1e[4]>13) {next;}
				 while(my $line=<RF>){
				 	 chomp($line);
				 	 if($.==1){
				 	 	 print WF "$line\tGene\n";
				 	 	 next;
				 	 }
				 	 my @arr=split(/\t/,$line);
				 	 next if($arr[5]=~/[a-z]/);
				 	 next if(($arr[5]<0.2) && ($arr[5]> -0.2));
				 	 foreach my $posline(@posline){
				 	 	 chomp($posLine);
				 	 	 my @posArr=split(/\t/,$posLine);
				 	 	 if($arr[1] eq $posArr[1]){
				 	 	 	  unless(($arr[2]>$posArr[3]) || ($arr[3]<$posArr[2])){
				 	 	 	  	print WF "$line\t$posArr[0]\n";
				 	 	 	  }
				 	 	 }
				 	 }
				 }
				 close(WF);
				 close(RF);
			}
		}
		close(SUB);
	}
}



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
