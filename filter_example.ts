import { FilterType, Filter, Butterworth } from "./butterworth_filter"
import * as fs from 'fs';

var fileDir: string = "E:\\developer\\SleepProject\\RingAlgorithm\\examples";
var fileName: string = "test1.txt";
var outputFileName: string = "jsoutput.txt";

var rawData: number[] = new Array();

export const Sleep = (ms: number) => {
    return new Promise(resolve => setTimeout(resolve, ms))
}

function readContent(dir: string, fName: string) {
    let fPath: string = dir.concat("\\").concat(fName);
    return fs.readFileSync(fPath, "utf8");
}

rawData = readContent(fileDir, fileName).split("\n").slice(0, -1).map(function (item) {
    return parseFloat(item.trim());
});
// console.log(rawData.length);

// resize to 2D array
var rawData2D: number[][] = new Array(rawData.length);
for (var i: number = 0; i < rawData.length; i++) {
    rawData2D[i] = new Array(1);
}
for(const index in rawData) {
    rawData2D[index][0] = rawData[index];
}
// console.log(rawData2D)

//filter configuration
var samplingRate: number = 200;
var filterOrder: number = 2;
var type: FilterType = FilterType.Bandpass;
var cutoff: number[] = [0.2, 0.4];

//create filter
var filt: Butterworth = new Butterworth(samplingRate, type, filterOrder, cutoff);
console.log(filt.Coefficients);

//apply filter
// var dataFilt: number[][] = Filter.filter(rawData, filt);

//apply filt filt
//filter滤波有明显的延迟，filtfilt滤波延时比filter小很多
//优先用filtfilt
var dataFiltFilt: number[][] = Filter.filtfilt(rawData2D, filt);
// console.log(dataFiltFilt);

var dataFiltFiltString = dataFiltFilt
    .map((item) => {
        return item;
    })
    .join("\n");

fs.writeFile(fileDir.concat("\\").concat(outputFileName), dataFiltFiltString, function (err) {
    if (err) {
        return console.error(err);
    }
    console.log("Created data file.");
});
