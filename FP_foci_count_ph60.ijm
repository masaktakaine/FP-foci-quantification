//出芽酵母細胞内のAde4-GFP fociを持つ細胞の割合をカウントするためのマクロ  2022.05.12
//位相差像から細胞の輪郭を抽出する
// Nikon Ti2-E@A103, 60倍位相差レンズを使用して位相差像とGFP/RFP蛍光像を取得することを想定
// textファイルをファイル数分作るとともに、視野中の全細胞数とfociを持つ細胞数を一つのファイルにまとめる
// FindMaxima関数のprominenceを入力するように変更

//個々のfociのxy座標をFindMaximaで取得した後、各座標でROIを作って平均輝度を計測
macro "Ade4-FP_foci_count" {

#@ String(label="Date of experiments, e.g., 2022-02-05") edate1
#@ int(label="Noise tolerance for FindMaxima", 3000) prom  // FindMaxima関数で使用する prominence = noise tolerance を入力
#@ File (label="Choose source Folder", style="directory") dirS0  // directory名の最後に"/"がつかない
#@ File (label="Choose destination Folder", style="directory") dirD0 // directory名の最後に"/"がつかない
#@ String(label="Hide/Show the active image? The Show slows the analysis.", choices={"hide","show"}, style="radioButtonHorizontal") sbm

setBatchMode(sbm); // hides the active image, required ImageJ 1.48h or later
dirS = dirS0 + File.separator; // File.separatorはOSに応じたファイル区切り文字に対応
dirD1 = dirD0 + File.separator;

//setOption("ExpandableArrays", true);  // expandableArrayオプションを使用すると配列のサイズを指定する必要が無い、ImageJ1.53以降は自動的にarrayがexpandするのでこのオプションは不要
edate = " "+edate1; // 日付の前にスペースを入れておくとexcelで勝手に変換されない

imagefilelist = getFileList(dirS);
file_name = newArray; 	//ファイル名を格納するための配列を作成
total_cell_number = newArray;	//各ファイル中の総細胞数を格納するための配列
foci_cell = newArray; 	//各ファイル中のfociを持つ細胞数を格納するための配列
pct_foci_cell = newArray; 	//各ファイル中のfociを持つ細胞の%を格納するための配列
prominence = newArray;
defocus = newArray;  // 解析に適さない画像をチェックするための配列、空にしておく

//exp_date = newArray;
dirD = dirD1 +  "/"+ "prom_" + prom;
File.makeDirectory(dirD);
dirBF = dirD + "/BF/";				//明視野画像のためのフォルダBFを作成
File.makeDirectory(dirBF);				//

dirGreen = dirD + "/green/";				//緑蛍光画像のためのフォルダGreenを作成
File.makeDirectory(dirGreen);				//

dirMer = dirD + "/merge/";				//マージ画像のためのフォルダmergeを作成
File.makeDirectory(dirMer);	

dirDR = dirD + "/Drawings/";				//Drawingsフォルダを作成してその中にマスク画像やROIを保存
File.makeDirectory(dirDR);			//

dirCSV = dirD +"/" + edate1 + "_prom" + prom + "_csv/";
File.makeDirectory(dirCSV);		// csv fileを保存するためのフォルダ

for (i = 0; i < imagefilelist.length; i++) {
   currFile = dirS + imagefilelist[i];
    if((endsWith(currFile, ".nd2"))||(endsWith(currFile, ".oib"))||(endsWith(currFile, ".zvi"))) { // process if files ending with .oib or .nd2, or .zvi
//   open(currFile); 								// Bio-Formatsのwindowを逐一表示、OKを押す必要がある
//  	run("Bio-Formats Windowless Importer", "open=currFile");
		run("Bio-Formats Macro Extensions"); 			//画像の明るさとコントラストはautoscaleされた状態、元々のピクセル値は不変
		Ext.openImagePlus(currFile)}
		else if ((endsWith(currFile, ".tif"))||(endsWith(currFile, ".tiff"))) {// process if files ending with .tif or .tiff
			open(currFile);   // 拡張子がtifのhyperstackファイルの場合は普通にopenする
		}

date = newArray;			//テーブルfoci_dataに日時を格納するための配列
image_file = newArray;  	//テーブルfoci_dataにファイル名を格納するための配列
cell_number =newArray;		//テーブルfoci_dataに細胞の番号を格納するための配列
obserbation = newArray;		//テーブルfoci_dataにobservationの通し番号を格納するための配列
foci_serial = newArray;		//テーブルfoci_dataにfociの通し番号を格納するための配列
foci_meanint =newArray;		//テーブルfoci_dataにfociの平均輝度を格納するための配列
foci_x = newArray;			//テーブルfoci_dataにfociのx座標を格納するための配列
foci_y = newArray;			//テーブルfoci_dataにfociのy座標を格納するための配列
cell_area = newArray;		//テーブルfoci_dataに細胞の面積を格納するための配列
prominences = newArray;		//テーブルfoci_dataにprominenceを格納するための配列
withfoci = newArray;		//テーブルfoci_dataにfociの有無をTRUE or FALSEで格納
defocuses = newArray;  // 解析に適さない画像をチェックするための配列、空にしておく

print("\\Clear"); 							// Log windowをリセット
title = getTitle();
title_s = replace(title, "\\.nd2", "");  // titleから拡張子.nd2を除く

imagewidth = getWidth();   // 画像の縦と横のサイズを取得
imageheight= getHeight();

run("Z Project...", "projection=[Max Intensity]");
title2 = getTitle();
run("Split Channels");

c1 = "C1-" + title2; // channel#1：Ade4-GFPのスタック画像
c2 = "C2-" + title2; // channel#2：60倍位相差明視野像

selectWindow(c2);
BFID = getImageID();
run("Duplicate...", " ");
run("Gaussian Blur...", "sigma=1");   //100倍では2
run("Subtract Background...", "rolling=25");  // 100倍では42

setAutoThreshold("Li dark");
setOption("BlackBackground", true);
run("Convert to Mask");

run("Dilate");
//run("Dilate");  //60倍ではdialate 1回

run("Skeletonize");
run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune_0");
run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune_0");  // ph100の場合、1回だけだと"Tagged skeleton"の画像が生成されないことがあるので2回繰り返す
selectWindow("Tagged skeleton");
run("Convert to Mask");
rename("debranch_" + title2);

run("Dilate"); //60倍ではdialate 1回
run("Close-");  //60倍ではclose1回

c3 = "debranch_" + title2;

run("Duplicate...", " ");
run("Fill Holes");
rename("filled_" + title2);

c4 = "filled_" + title2;

imageCalculator("Subtract create", c4,c3);
selectWindow("Result of " + c4);
rename("segmented_" + title2);
run("Set Measurements...", " area redirect=None decimal=3");
run("Analyze Particles...", "size=400-3000 pixel show=Outlines display exclude clear add"); // ph100では667-5000, 細胞の輪郭を検出してROI managerに送る（add）

//waitForUser("Select the image for spot counting"); どの画像を使用するか選択させる場合
sno = 0; // 観察数を数えるための変数： Serial No. of observation, 最終的に生成されるfoci_dataの行数に等しい、総foci数+foci持たない細胞数
cwf = 0; // Fociを持つ細胞を数えるための変数: Cell With Foci
snf = 0; // 各ファイルにおける総foci数を数えるための変数: Serial No. of Foci

selectWindow(c1);							// 検出した細胞の輪郭（ROI）をAde4-mNG imageに重ねる
roiManager("Set Color", "white");
roiManager("Show None");
roiManager("Show All");
roiManager("Save", dirDR + "ROIs-" + title2 + ".zip"); 
run("Clear Results");

//run("Set Measurements...", "mean redirect=None decimal=3");	// multi-point selectionで検出したfociのmean intensityを測定するため
Table.create("cell_info");										//各細胞のfociの有無、面積を書き込むためのテーブルを作っておく

for(k=0; k<roiManager("count"); k++) {   // ROI managerでROIを一つずつactiveにして解析
	run("Clear Results");
	roiManager("select", k);
	Roi.setStrokeColor(255*random,255*random,255*random);
	area = getValue("Area");		//選択したROIの面積を計算して変数areaに格納		
	run("Find Maxima...", "prominence=" + prom + " exclude output=Count");// prominence = noise tolerance, excludeはedgeにおけるmaximaを除く,maximaの数をカウント
	
	nf_k = getResult("Count", 0);	//0行目に表示されている、maximaの数=fociの数を変数nf_kに格納、Tableを一つしか使っていない（Resultsのみ）場合はTable.getでもgetResultsでも同じ結果になる
									// number of foci in k-th cell
	run("Clear Results");			//fociがある場合はMeasureするのでResultsテーブルをリセットする
	if(nf_k >= 1){						// k番目の細胞にfociが存在する場合
		//selectWindow("cell_info");						// テーブルcell_infoを選択
		//Table.set("Cell with foci", k, 1);				//その細胞がfociを持つ場合、cell with fociに1を入れる
		cwf = cwf + 1;									// fociを持つならばcwfを一つ増やす
		
		run("Find Maxima...", "prominence=" + prom + " exclude output=List");	//maximaを検出して、それらのxy座標をResultsテーブルに表示
		for (l=0; l< nf_k; l++){						// k番目の細胞におけるnf_k個のfociについての解析開始
			obserbation[sno] = sno+1;
			foci_serial[sno] = snf+1;					// fociの通し番号の計算
			cell_number[sno] = k+1;						// snf番目のfociがk番目の細胞にあることを記述
			image_file[sno] = title_s;					//拡張子を除いたファイル名を格納
			
			foci_x[sno] = getResult("X", l);			// snf番目のfociのX座標をResultsテーブルから取得
			foci_y[sno] = getResult("Y", l);			// snf番目のfociのY座標をResultsテーブルから取得
			cell_area[sno] = area;
			prominences[sno] = prom;
			date[sno] = edate;
			withfoci[sno]="TRUE";

			//makeOval(foci_x[snf]-2,foci_y[snf]-2,4,4);	// snf番目のfociを囲むような直径4pixelの円のROIを作る、指定するxy座標は円を囲む四角形の左上の座標
														// ph60では直径2pixel: makeOval(foci_x[snf]-1,foci_y[snf]-1,2,2)
			run("Specify...", "width=2 height=2 x="+foci_x[sno]+" y="+foci_y[sno]+" slice=1 oval centered"); //これでも良いかもしれない、centeredを指定すれば座標をそのまま使用できそう
			foci_meanint[sno] = getValue("Mean");		// 円内の平均輝度を取得してfoci_meanintに格納
	
			Table.update;
			sno = sno+1;
			snf= snf+1;								// snfを一つ増やす
		}
	}else {												// nf_k = 0,すなわちk番目の細胞にfociが存在しない場合
		//selectWindow("cell_info");						// テーブルcell_infoを選択
		//Table.set("Cell with foci", k, 0);				//その細胞がfociを持たない場合、cell with fociに0を入れる
		//Table.set("Area", k, area);						// nf_kの値に関わらず、テーブルcell_infoにk番目の細胞の面積を書き込む
		
		obserbation[sno] = sno+1;
		foci_serial[sno] = "NaN"; 					// fociが無い場合
		cell_number[sno] = k+1;
		image_file[sno] = title_s;					//拡張子を除いたファイル名を格納
		foci_x[sno] = "NaN";			// 検出されなかったことを意味するNaN
		foci_y[sno] = "NaN";			// 検出されなかったことを意味するNaN
		cell_area[sno] = area;
		foci_meanint[sno] = "NaN";		// 検出されなかったことを意味するNaN
		prominences[sno] = prom;
		date[sno] = edate;
		withfoci[sno]="FALSE";
		Table.update;
		sno=sno+1;
}
}
//Table.set("# of cells", k, k);							//テーブルcell_infoの最後の行に総細胞数kを書き込む,forループを抜けた後なので、kは総細胞になっていることに注意
//Table.set("# of cells w/ foci", k, cwf);				//テーブルcell_infoの最後の行にfociを持っている細胞数を書き込む
//Table.set("Filename",k,title_s);							
//Table.set("Prominence", k, prom);


file_name[i] = title_s;									// i番目のファイルの名前を格納, 拡張子無し
total_cell_number[i] = k;								// i番目のファイルにおける総細胞数を格納
foci_cell[i] = cwf;										// i番目のファイルにおけるfociを持つ細胞数を格納
pct_foci_cell[i] = 100*cwf/k;							// i番目のファイルにおけるfociを持つ細胞の割合（%）を格納
prominence[i] = prom;

saveAs("Tiff", dirGreen+getTitle());

c1t = replace(c1, "\\.nd2", "\\.tif");					// C1の方は拡張子がtifになっているので、それに合わせた名前を作成

run("Merge Channels...", "c4=["+c2+"] c6=["+c1t+"] keep ignore"); // run()のoptionにマクロ中の変数を入力する場合は、""と++で囲む

selectWindow("RGB");
rename("BFM-" + title);
RGID = getImageID();
saveAs("Tiff", dirMer+getTitle()); // merge画像を保存
close();

selectWindow("cell_info");
//saveAs("Results", dirTX + title_s + ".txt");
run("Close");	// 画像以外のwindowを閉じる時に使用

selectImage(BFID);
saveAs("Tiff", dirBF+getTitle());
selectWindow("Drawing of segmented_"+ title2);
saveAs("Tiff", dirDR+getTitle());
close();

// fociの平均輝度を測定したROIを描画する
newImage(title_s+"_foci_location", "8-bit", imagewidth, imageheight, 1);
setColor(0,0,0);
fill();
//setLineWidth(1);
for (m = 0; m < foci_x.length; m++) {
	setColor(255,255,255);
	fillOval(foci_x[m]-2, foci_y[m]-2, 4, 4);
}
saveAs("Tiff", dirDR+getTitle());
close();

Array.show("foci_data", date,image_file,cell_number, obserbation,withfoci, foci_serial, foci_meanint, foci_x, foci_y, cell_area, prominences, defocuses); //配列を独立したwindowに表示する
//Array.show("foci_data(row numbers)", image_file, cell_number, foci_serial, foci_meanint, foci_x, foci_y); //配列を独立したwindowに表示する、タイトルが(row numbers)で終わると最初の列が1から始まる行のインデックスになる
selectWindow("foci_data");					// 個々のfociの情報について記述した表になる
saveAs("Results", dirCSV + edate1 + "_"+ title_s + "_prom" + prom + ".csv"); //以前は拡張子xlsにしていたが、紛らわしいのでcsvにする 
run("Close");	// 画像以外のwindowを閉じる時に使用
run("Close All");
//close("*.xls");
//close("*.txt");

}
    
  //Array.show("foci_stat(row numbers)", image_file, cell_number, foci_cell, pct_foci_cell, prominence, exp_date); //配列を独立したwindowに表示する、タイトルが(row numbers)で終わると最初の列が1から始まる行のインデックスになる
 Array.show("foci_stat(row numbers)", file_name, total_cell_number, foci_cell, pct_foci_cell, prominence, defocus); //配列を独立したwindowに表示する、タイトルが(row numbers)で終わると最初の列が1から始まる行のインデックスになる
																	
    selectWindow("foci_stat");
    saveAs("Results", dirD +"/"+ edate1 + "_prom" + prom + "_foci_stat.csv"); 
   
	run("Close");	// 画像以外のwindowを閉じる時に使用
    run("Clear Results");						// Results windwoをリセット
	print("\\Clear"); 							// Log windowをリセット
	roiManager("reset");						// ROI managerをリセット

  showMessage(" ", "<html>"+"<font size=+2>Process completed<br>");
}