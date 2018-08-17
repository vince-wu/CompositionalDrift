numUniqueMonomers = 2;
//================================================================================================================================
//SIMULATION AND DYNAMIC VISUALIZATION
//================================================================================================================================
function simulate() {
	colorArray = ["#FF6600", "#FCD202", "#B0DE09", "#0D8ECF", "#2A0CD0", 
	"#CD0D74", "#CC0000", "#00CC00", "#0000CC", "#DDDDDD", "#999999", "#333333", "#990000"];
	//Get all relevant inputs
	numUniqueMonomers = parseInt(document.getElementById("numUniqueMonomers").value);
	var totalNumMonomers = parseInt(document.getElementById("totalNumMonomers").value);
	var mRatio = parseInt(document.getElementById("mRatio").value);
	var conversion = parseInt(document.getElementById("conversion").value);
	var numRowsToShow = parseInt(document.getElementById("numRowsToShow").value);
	var graphTypeObj = document.getElementById("graph1Type");
	var graphType = graphTypeObj.options[graphTypeObj.selectedIndex].value;
	//Initial variable calculations
	//console.log("Type of monomer1Ratio: ", typeof(monomer1Ratio));
	polymerLength = Math.floor(mRatio * conversion / 100);
	var numPolymers = Math.floor(totalNumMonomers / mRatio);
	var monomerRatioList = getMonomerRatios(numUniqueMonomers);
	console.log("monomerRatioList: ", monomerRatioList);

	//Set up visualizing canvas
	var canvas = document.getElementById("visual");
	var ctx = canvas.getContext("2d");
	var winWidth = window.innerWidth;
	var canvasWidth = parseInt(winWidth * 0.95);
	var squareLength = Math.min((canvasWidth - 10) / polymerLength, 15);
	canvas.setAttribute("width", canvasWidth);
	canvas.setAttribute("height", squareLength * numRowsToShow + 10);
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	ctx.beginPath();
	ctx.translate(0.5, 0.5);
	ctx.lineWidth="1";
	ctx.strokeStyle = "black";
	//var rrList = [[monomer1RR, 1], [1, monomer2RR]];
	var monomerAmountsList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	initialMonomerAmountList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	console.log("initialMonomerAmountList1: ", initialMonomerAmountList);
	rrList = getReactivityRatios(numUniqueMonomers);
	//console.log(monomerAmountsList);
	//Initiate chains, all polymer chains are represented as arrays and stored in another array, polymerArray
	polymerArray = [];
	var currNumPolymers;
	for (currNumPolymers = 0; currNumPolymers < numPolymers; currNumPolymers++) {
		var initChoices = [];
		var initWeightList = [];
		//if it is a 2-monomer system, use the Mayo Lewis Equation to Determine Starting Monomers
		//Weights that determine probabilty of monomer starting the chain are calculated
		if (numUniqueMonomers == 2) {
			var f1 = monomerAmountsList[0];
			var f2 = monomerAmountsList[1];
			var r1 = rrList[0][0];
			var r2 = rrList[1][1];
			var initWeight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2);
			initChoices.push(1);
			initChoices.push(2);
			initWeightList.push(initWeight);
			initWeightList.push(1 - initWeight);
		} else {
		//If it is a 3-monomer or more system, use initial starting ratios to determing starting monomers
			var initChoices = [];
			var initWeightList = [];
			for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
				initChoices.push(monomerID);
				var initWeight = monomerAmountsList[monomerID - 1];
				initWeightList.push(initWeight);	
			}
		}
		//Use a weighted random selector to choose the starting monomer
		var startingMonomer = parseInt(weightedRand(initChoices, initWeightList));
		//console.log("typeof startingMonomer: ", typeof(startingMonomer));
		//console.log("startingMonomer: ", startingMonomer);
		//Create a new polymer chain with the selected starting monomer
		polymerArray.push([startingMonomer]);
		//Remove one of that monomer from the pool
		monomerAmountsList[startingMonomer - 1] --;
		//Visualize the monomer being added
		// if (currNumPolymers < numRowsToShow) {
		// 	color = colorArray[startingMonomer - 1];
		// 	ctx.beginPath()
		// 	ctx.fillStyle = color;
		// 	ctx.rect(0, currNumPolymers*squareLength, squareLength, squareLength);
		// 	ctx.fill();
		// 	ctx.stroke();
		// 	sleep(100);
		// }
	}
	//console.log("polymerArray: ", polymerArray);
	//Propogate chains
	//First loop: propagate until the polymer length reaches expected value
	for (var currLength = 1; currLength < polymerLength; currLength++) {
		//Second loop: during each cycle of propagation, add one monomer to each of the polymer chains
		for (var i = 0; i < polymerArray.length; i++) {
			//console.log("polymerArray: ", polymerArray);
			var polymer = polymerArray[i];
			var propagationChoices = [];
			var propagationWeightList = [];
			//Third loop: cycle through each unique monomer and calculate probability weight of the monomer to be added to 
			// the growing polymer chain
			for (var monomerID = 1; monomerID < numUniqueMonomers + 1; monomerID++) {
				//console.log("rrList: ", rrList);
				//console.log("polymer: ", polymer);
				//console.log("last monomer: ", polymer[polymer.length - 1] - 1);
				var reactivityRatio = rrList[polymer[polymer.length - 1] - 1][monomerID - 1];
				var propWeight = monomerAmountsList[monomerID - 1] * reactivityRatio;
				//console.log("monomerID: ",monomerID, "propWeight: ", propWeight);
				propagationChoices.push(monomerID);
				propagationWeightList.push(propWeight);
			}
			//Using weights, randomly select propagating monomer and add it to the polymer chain array
			//console.log("propagationSpec:", propagationSpec);
			var nextMonomer = parseInt(weightedRand(propagationChoices, propagationWeightList));
			polymer.push(nextMonomer);
			//Remove that monomer from the reactant pool
			monomerAmountsList[nextMonomer - 1]--;
			// if (i < numRowsToShow) {
			// 	color = colorArray[nextMonomer - 1];
			// 	length = squareLength;
			// 	console.log("currLength: ", currLength);
			// 	xfactor = currLength;
			// 	yfactor = i;
			// 	function draw() {
			// 		console.log("currLength2: ", currLength)
			// 		ctx.beginPath()
			// 		ctx.fillStyle = color;
			// 		ctx.rect(xfactor*length, yfactor*length, length, length);
			// 		ctx.fill();
			// 		ctx.stroke();
			// 	}
			// 	requestAnimationFrame(draw)
			// }
		}
	}
	visualize(canvas, ctx, polymerArray, polymerLength, numRowsToShow, squareLength);
	setGraph(graphType);
};

//================================================================================================================================
//FUNCTIONS FOR RETRIEVING USER INPUT DATA
//=================================================================================================================================

//Parses HTML document, returns list of monomer ratios
function getMonomerRatios(numUniqueMonomers) {
	var monomerRatioList = [];
	var inputTable = document.getElementById("inputTable");
	for (var index = 0; index < numUniqueMonomers; index++) {
		var row = inputTable.rows[index];
		var cell = row.cells[5];
		var inputObj = cell.children[0];
		var monomerRatio = parseFloat(inputObj.value); 
		monomerRatioList.push(monomerRatio);
	}
	return monomerRatioList;
};

//Parses HTML document, returns a list of lists of reactivity ratios

function getReactivityRatios(numUniqueMonomers) {
	var rrList = [];
	var inputTable = document.getElementById("inputTable");
	if (numUniqueMonomers == 2) {
		var row1 = inputTable.rows[0];
		var cell1 = row1.cells[7];
		var inputObj1 = cell1.children[0];
		var rr1 = parseFloat(inputObj1.value);
		var singleRRList1 = [rr1, 1];
		var row2 = inputTable.rows[1];
		var cell2 = row2.cells[7];
		var inputObj2 = cell2.children[0];
		var rr2 = parseFloat(inputObj2.value);
		var singleRRList2 = [1, rr2];
		rrList = [singleRRList1, singleRRList2];
		return rrList;
	}
	for (var index = 0; index < numUniqueMonomers; index++) {
		var singleRRList = []; 
		rrList.push(singleRRList);
		for (var index2 = 0; index2 < numUniqueMonomers; index2++) {
			var row = inputTable.rows[index2];
			var cell = row.cells[7 + 2 * index];
			var inputObj = cell.children[0];
			var reactivityRatio = parseFloat(inputObj.value);
			singleRRList.push(reactivityRatio);
		}
	}
	return rrList;
};

//==================================================================================================================================
//CREATING ITERATIVE INPUTS
//==================================================================================================================================

function createInputs(inputNum) {
	inputNum = parseFloat(inputNum);
	if (inputNum < 2) {
		return;
	}
	var inputTable = document.getElementById("inputTable");

	//Change HTML Objects for Monomer Ratios

	//Delete all current objects

	for (var index = 0; index < numUniqueMonomers; index++) {
		var row = inputTable.rows[index];
		row.deleteCell(4);
		row.deleteCell(4);
	}

	//Iteratively create objects

	for (var index = 0; index < inputNum; index++) {
		var labelElement = document.createElement("LABEL");;
		labelElement.innerHTML = "Monomer " + (index + 1) + " Ratio:";
		var inputElement = document.createElement("INPUT");
		inputElement.setAttribute("type", "number");
		inputElement.setAttribute("value", 1);
		var row = inputTable.rows[index];
		//console.log("row: ", "row"+index);
		var labelCell = row.insertCell(4);
		labelCell.appendChild(labelElement);
		var inputCell = row.insertCell(5);
		inputCell.appendChild(inputElement);
	}

	//Change HTML Objects for Monomer Reactivities

	//Delete all objects

	for (var index = 0; index < numUniqueMonomers; index++) {
			var row = inputTable.rows[index];
			if (numUniqueMonomers == 2) {
				row.deleteCell(-1);
				row.deleteCell(-1);
			} else {
				for (var index2 = 1; index2 <= numUniqueMonomers; index2++) { 
					row.deleteCell(-1);
					row.deleteCell(-1);
				}
			}
	}

	//Iteratively create objects 	
	if (inputNum > 2) {
		for (var index = 0; index < inputNum; index++) {
			var row = inputTable.rows[index];
			for (var index2 = 1; index2 <= inputNum; index2++) {
				var labelElement = document.createElement("LABEL");
				var labelText = index2 + "-" + (index + 1) + " Reactivity";
				labelElement.innerHTML = labelText;
				var inputElement = document.createElement("INPUT");
				inputElement.setAttribute("type", "number");
				inputElement.setAttribute("value", 1);
				var labelCell = row.insertCell(-1);
				labelCell.appendChild(labelElement);
				var inputCell = row.insertCell(-1);
				inputCell.appendChild(inputElement);
			}
		}
	}
	if (inputNum == 2) {
		for (var index = 0; index < 2; index++) {
			var row = inputTable.rows[index];
			var labelElement = document.createElement("LABEL");
			var labelText = "Monomer " + (index + 1) + " Reactivity";
			labelElement.innerHTML = labelText;
			var inputElement = document.createElement("INPUT");
			inputElement.setAttribute("type", "number");
			inputElement.setAttribute("value", 1);
			var labelCell = row.insertCell(-1);
			labelCell.appendChild(labelElement);
			var inputCell = row.insertCell(-1);
			inputCell.appendChild(inputElement);
		}
	}

	//Update numUniquemonomers

	numUniqueMonomers = inputNum;

};

//=========================================================================================================================================
//POLYMER ANALYSIS FUNCTIONS
//=========================================================================================================================================

function getMonomerAmounts(monomerRatioList, totalNumMonomers) {
	var monomerAmountsList = [];
	var totalWeight = sumArray(monomerRatioList);
	for (var i = 0; i < monomerRatioList.length; i++) {
		var monomerRatio = monomerRatioList[i];
		var weight = monomerRatio / totalWeight;
		var monomerAmount = Math.ceil(totalNumMonomers * weight);
		//console.log("monomerRatio: ", monomerRatio);
		monomerAmountsList.push(monomerAmount);
	}
	return monomerAmountsList;
};
//Returns chartData formatted information for percent composition of each unique monomer at each index.
function getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength) {
	var fullCompositionList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singleCompositionList = [];
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			var monomerCount = 0;
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					monomerCount++;
				}
			}
			var monomerRatio = monomerCount / polymerArray.length;
			singleCompositionList.push(monomerRatio);
		}
		var xyDataPair = [createRangeArray(polymerLength), singleCompositionList];
		fullCompositionList.push(xyDataPair);
	}
	chartData = convertToChartData(fullCompositionList, numUniqueMonomers, polymerLength);
	return chartData;
};
function getPercentageMonomer(polymerArray, numUniqueMonomers, polymerLength, initialMonomerAmountList) {
	var fullPercentageList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singlePercentageList = [];
		var monomerRemaining = initialMonomerAmountList[monomerID -1];
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					monomerRemaining--;
				}
			}
			//console.log("initialMonomerAmountList: ", initialMonomerAmountList);
			var monomerPercentage = monomerRemaining / initialMonomerAmountList[monomerID -1];
			singlePercentageList.push(monomerPercentage);
		}
		var xyDataPair = [createRangeArray(polymerLength), singlePercentageList];
		fullPercentageList.push(xyDataPair);
	}
	chartData = convertToChartData(fullPercentageList, numUniqueMonomers, polymerLength);
	return chartData;
};
function getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, monomerID){
	//console.log("monomerID: ", monomerID);
	//console.log("typeof monomerID: ", typeof(monomerID));
	var fullSeparationData = new Array(polymerLength + 1).fill(0)
	var largestBlock = 0;
	for (var i = 0; i < polymerArray.length; i++) {
		polymer = polymerArray[i];
		var numConsecutive = 0;
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			//console.log("numConsecutive: ", numConsecutive);
			currMonomerID = polymer[positionIndex];
			//console.log("currMonomerID: ", currMonomerID);
			if (currMonomerID != monomerID && numConsecutive > 0) {
				fullSeparationData[numConsecutive] += numConsecutive;
				if (numConsecutive > largestBlock) {
					largestBlock = numConsecutive;
				}
				numConsecutive = 0;
			}
			if (currMonomerID == monomerID) {
				numConsecutive += 1;
			}
			if (positionIndex == polymerLength - 1) {
				if (numConsecutive != 0) {
					fullSeparationData[numConsecutive] += numConsecutive;
				}
			}
		}
	}
	console.log("largestBlock: ", largestBlock);
	fullSeparationData = fullSeparationData.slice(0, largestBlock + 2);
	var chartData = [];
	for (var blockSize = 0; blockSize < fullSeparationData.length; blockSize++) {
		var obj = {};
		var xID = "blockSize" + monomerID;
		var xCounts = "counts" + monomerID;
		obj[xID] = blockSize;
		obj[xCounts] = fullSeparationData[blockSize];
		chartData.push(obj);
	}
	return chartData;		
};
function getPolymerCompostion(polymerArray, numUniqueMonomers, polymerLength) {
	var fullCompositionList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singleCompositionList = [];
		var totalMonomerCount = 0;
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					totalMonomerCount++;
				}
			}
			var polymerComposition = totalMonomerCount / (polymerArray.length * (positionIndex + 1));
			singleCompositionList.push(polymerComposition);
		}
		var xyDataPair = [createRangeArray(polymerLength), singleCompositionList];
		fullCompositionList.push(xyDataPair);
	}
	chartData = convertToChartData(fullCompositionList, numUniqueMonomers, polymerLength);
	return chartData;
};
function convertToChartData(fullList, numUniqueMonomers, polymerLength) {
	chartData = [];
	for (var i = 0; i <= polymerLength; i++) {
		var obj = {};
		for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
			var xID = "x" + monomerID;
			var yID = "y" + monomerID;
			//console.log("xID: ", xID);
			obj[xID] = fullList[monomerID - 1][0][i];
			obj[yID] = fullList[monomerID - 1][1][i];
		}
		chartData.push(obj);
	}
	return chartData;
};

//==========================================================================================================================================
//GRAPHING FUNCTIONS
//==========================================================================================================================================

function showDemoGraph() {
	colorArray = ["#FF6600", "#FCD202", "#B0DE09", "#0D8ECF", "#2A0CD0", 
	"#CD0D74", "#CC0000", "#00CC00", "#0000CC", "#DDDDDD", "#999999", "#333333", "#990000"];
	// XY CHART
	chart = new AmCharts.AmXYChart();
	chart.dataProvider = demoData;
	chart.startDuration = 0;
	chart["export"] = {"enabled": true};

	//X
	var xAxis = new AmCharts.ValueAxis();
	xAxis.title = "Monomer Index";
	xAxis.position = "bottom";
	xAxis.dashLength = 1;
	xAxis.axisAlpha = 0;
	xAxis.autoGridCount = true;
	xAxis.minimum = 0;
	xAxis.maximum = 100;
	chart.addValueAxis(xAxis);

	// Y
	var yAxis = new AmCharts.ValueAxis();
	yAxis.position = "left";
	yAxis.title = "Monomer Occurence";
	yAxis.dashLength = 1;
	yAxis.axisAlpha = 0;
	yAxis.autoGridCount = true;
	yAxis.minimum = 0;
	yAxis.maximum = 1;
	chart.addValueAxis(yAxis);

	//Remove all old graphs
	var graphList = chart.graphs;
	for (var i = 0; i < graphList.length; i++) {
		chart.removeGraph(graphList[i]);
	}
	//Add graphs depending on number of unique monomers
	for (var graphNumber = 1; graphNumber <= 2; graphNumber++) {
		var graph = new AmCharts.AmGraph();
	  	graph.lineColor = colorArray[graphNumber - 1];
		graph.balloonText = "x:[[x]] y:[[y]]";
		graph.xField = "x" + graphNumber;
		graph.yField = "y" + graphNumber;
		graph.lineAlpha = 0;
		graph.bullet = "triangleDown";
		graph.title = "Monomer" + graphNumber;
		chart.addGraph(graph);
	}
	//LEGEND
	var legend = new AmCharts.AmLegend();
	legend.position = "right";
	legend.align = "center";
	legend.markerType = "square";
	legend.useGraphSettings = true;
	chart.addLegend(legend);

	// CURSOR
	var chartCursor = new AmCharts.ChartCursor();
	chart.addChartCursor(chartCursor);

	// SCROLLBAR

	var chartScrollbar = new AmCharts.ChartScrollbar();
	chartScrollbar.scrollbarHeight = 5;
	chartScrollbar.offset = 15
	chart.addChartScrollbar(chartScrollbar);
	// WRITE
	chart.write("chartdiv");
};

function createHistChart(chartData, monomerID) {
	chart2 = new AmCharts.AmSerialChart();
	var obj = chartData[0];
	keyList = Object.keys(obj);
	console.log("keyList: ", keyList);
	//Histogram
	chart2.dataProvider = chartData;
	chart2.categoryField = keyList[0];
	chart2.startDuration = 1;
	chart2.sequencedAnimation = false; 

	//Value Axis
	var valueAxis = new AmCharts.ValueAxis();
	valueAxis.position = "bottom";
	valueAxis.minimum = 0;
	valueAxis.title = "Normalized Separation";
	chart2.addValueAxis(valueAxis);

	//X Axis

	chart2.categoryAxis.title = "Monomer " + monomerID + " Block Size";

	//LEGEND

	var legend = new AmCharts.AmLegend();
	legend.position = "bottom";
	legend.align = "center";
	legend.markerType = "square";
	legend.useGraphSettings = true;
	//chart2.addLegend(legend);

	//Graph
	var graph = new AmCharts.AmGraph();
	graph.type = "column";
	graph.fillAlphas = 1;
	graph.valueField = keyList[1];
	graph.balloonText = "[[value]]";
	chart2.addGraph(graph);
	chart2.write("chartdiv2");
};

function formatXYGraphs() {
	var chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);

	//Remove all old graphs
	var graphList = chart.graphs;
	//console.log("graphList.length", graphList.length);
	for (var i = 0; i < graphList.length; i++) {
		chart.removeGraph(graphList[i]);
	}
	console.log("lengthafter: ", graphList.length);
	//Add graphs depending on number of unique monomers
	for (var graphNumber = 1; graphNumber <= numUniqueMonomers; graphNumber++) {
		var graph = new AmCharts.AmGraph();
	  	graph.lineColor = colorArray[graphNumber - 1];
		graph.balloonText = "x:[[x]] y:[[y]]";
		graph.xField = "x" + graphNumber;
		graph.yField = "y" + graphNumber;
		graph.lineAlpha = 0;
		graph.bullet = "triangleDown";
		chart.addGraph(graph);
	}
	chart.removeGraph(graphList[0]);
};
function formatHistGraph(monomerID) {
	if (monomerID > numUniqueMonomers) {
		return;
	}
	monomerID = parseInt(monomerID);
	console.log("monomerID: ", monomerID);
	var chartData = getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, monomerID);
	var obj = chartData[0];
	keyList = Object.keys(obj);
	chart2.categoryField = keyList[0];
	chart2.graphs[0].valueField = keyList[1];
	chart2.dataProvider = chartData;
	chart2.categoryAxis.title = "Monomer " + monomerID + " Block Size";
	chart2.graphs[0].lineColor = colorArray[monomerID - 1];
	chart2.graphs[0].fillColors = colorArray[monomerID - 1];
	chart2.validateData();
	chart2.validateNow();
};
function setGraph(type) {
	//console.log("type: ", type);
	formatXYGraphs();
	var chartData;
	switch (type) {
		case "Monomer Occurences":
			chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);
			//console.log("chartData: ", chartData);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer Index";
			chart.valueAxes[1].title = "Monomer Occurence";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Percentage Monomer":
			chartData = getPercentageMonomer(polymerArray, numUniqueMonomers, polymerLength, initialMonomerAmountList);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer Index";
			chart.valueAxes[1].title = "Unreacted Monomer Remaining";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Polymer Compositions":
			chartData = getPolymerCompostion(polymerArray, numUniqueMonomers, polymerLength);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer index";
			chart.valueAxes[1].title = "Instantaneous Composition of Polymer Chain";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Monomer Separation":
			var histMonomer = parseInt(document.getElementById("hist1Monomer").value);
			chartData = getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, histMonomer);
			createHistChart(chartData, histMonomer);
			document.getElementById("chartdiv").style.display = "none";
			document.getElementById("chartdiv2").style.display = "block";
			chart2.invalidateSize();
			chart2.animateAgain();
			return;
	}
	document.getElementById("chartdiv2").style.display = "none";
	document.getElementById("chartdiv").style.display = "block";
	//console.log("got here");
	chart.invalidateSize();
};

//=================================================================================================================================
//VISUALIZING FUNCTIONS
//=================================================================================================================================

function visualize(canvas, ctx, polymerArray, polymerLength, numRowsToShow, length) {
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	ctx.beginPath();
	var blockList = [];
	for (var len = 0; len < polymerLength; len++) {
		console.log("got here!");
		// polymerIndex = 0;
		// var interval = setInterval(loop, 20);
		// function loop() {
		// 	console.log("polymerIndex: ", polymerIndex);
		// 	console.log("numRowsToShow: ", numRowsToShow);
		// 	console.log("len: ", len);
		// 	if (polymerIndex >= numRowsToShow) {
		// 		clearInterval(interval);
		// 		return;
		// 	}
		// 	//console.log("polymerArray: " ,polymerArray)
		// 	monomerID = polymerArray[polymerIndex][len];
		// 	var color = colorArray[monomerID - 1];
		// 	ctx.beginPath()
		// 	ctx.fillStyle = color;
		// 	ctx.rect(len*length, polymerIndex*length, length, length);
		// 	ctx.fill();
		// 	ctx.stroke();
		// 	polymerIndex += 1;
		// }
		for (var polymerIndex = 0; polymerIndex < numRowsToShow; polymerIndex ++) {
			monomerID = polymerArray[polymerIndex][len];
			var color = colorArray[monomerID - 1];
			var block = {
				bcolor: color,
				ulX: len * length,
				ulY: polymerIndex * length
			}
			blockList.push(block);
			//draw();
			//sleep(20)
			// function draw() {
			// 	ctx.beginPath()
			// 	ctx.fillStyle = color;
			// 	ctx.rect(len*length, polymerIndex*length, length, length);
			// 	ctx.fill();
			// 	ctx.stroke();
			// 	requestAnimationFrame(draw);
			// }
		}
	}
	var lastIndex = 0;
	function draw() {
		for (var index = 0; index < lastIndex; index++) {
			block = blockList[index];
			ctx.beginPath()
			ctx.fillStyle = block.bcolor;
			ctx.rect(block.ulX, block.ulY, length, length);
			ctx.fill();
			ctx.stroke();
		}
		if (lastIndex == blockList.length) {
			animation.pause();
		}
		if (lastIndex < blockList.length) {
			lastIndex += 1;
		}
		requestAnimationFrame(draw)
	}
	draw();
	// var ctx = canvas.getContext("2d");
	// ctx.fillRect(0,0,50,50);
	// ctx.fillRect(100,0,50,50);
	// ctx.fillRect(200,0,50,50);
}

function setCanvasScalingFactor() {
	return window.devicePixelRatio || 1;
}
//================================================================================================================================
//UTILITY FUNCTIONS
//================================================================================================================================
//Creates an array from [1, 2, 3,......, N]
function createRangeArray(length) {
	var array = [];
	for (var i = 1; i <= length; i++) {
		array.push(i);
	}
	return array;
};
//Sums the numerical contents of an array
function sumArray(array) {
	return array.reduce(sum, 0)
};
function sum(a, b) {
	return a + b;
};
function rand(min, max) {
    return Math.random() * (max - min) + min;
};
 
function weightedRand(list, weight) {
    var total_weight = weight.reduce(function (prev, cur, i, arr) {
        return prev + cur;
    });
     
    var random_num = rand(0, total_weight);
    var weight_sum = 0;
    //console.log(random_num)
     
    for (var i = 0; i < list.length; i++) {
        weight_sum += weight[i];
        weight_sum = +weight_sum.toFixed(2);
         
        if (random_num <= weight_sum) {
            return list[i];
        }
    }
    // end of function
};
//sleep function
function sleep(milliseconds) {
  var start = new Date().getTime();
  for (var i = 0; i < 1e7; i++) {
    if ((new Date().getTime() - start) > milliseconds){
      break;
    }
  }
}
//========================================================================================================================================
//DATA
//========================================================================================================================================

demoData = 
[{x1: 1, y1: 0.389, x2: 1, y2: 0.611},

{x1: 2, y1: 0.365, x2: 2, y2: 0.635},

{x1: 3, y1: 0.3755, x2: 3, y2: 0.6245},

{x1: 4, y1: 0.3905, x2: 4, y2: 0.6095},

{x1: 5, y1: 0.3635, x2: 5, y2: 0.6365},

{x1: 6, y1: 0.3785, x2: 6, y2: 0.6215},

{x1: 7, y1: 0.386, x2: 7, y2: 0.614},

{x1: 8, y1: 0.371, x2: 8, y2: 0.629},

{x1: 9, y1: 0.3945, x2: 9, y2: 0.6055},

{x1: 10, y1: 0.3705, x2: 10, y2: 0.6295},

{x1: 11, y1: 0.3985, x2: 11, y2: 0.6015},

{x1: 12, y1: 0.4015, x2: 12, y2: 0.5985},

{x1: 13, y1: 0.372, x2: 13, y2: 0.628},

{x1: 14, y1: 0.411, x2: 14, y2: 0.589},

{x1: 15, y1: 0.3755, x2: 15, y2: 0.6245},

{x1: 16, y1: 0.39, x2: 16, y2: 0.61},

{x1: 17, y1: 0.393, x2: 17, y2: 0.607},

{x1: 18, y1: 0.393, x2: 18, y2: 0.607},

{x1: 19, y1: 0.3915, x2: 19, y2: 0.6085},

{x1: 20, y1: 0.386, x2: 20, y2: 0.614},

{x1: 21, y1: 0.3915, x2: 21, y2: 0.6085},

{x1: 22, y1: 0.4015, x2: 22, y2: 0.5985},

{x1: 23, y1: 0.3855, x2: 23, y2: 0.6145},

{x1: 24, y1: 0.419, x2: 24, y2: 0.581},

{x1: 25, y1: 0.385, x2: 25, y2: 0.615},

{x1: 26, y1: 0.4125, x2: 26, y2: 0.5875},

{x1: 27, y1: 0.389, x2: 27, y2: 0.611},

{x1: 28, y1: 0.399, x2: 28, y2: 0.601},

{x1: 29, y1: 0.397, x2: 29, y2: 0.603},

{x1: 30, y1: 0.415, x2: 30, y2: 0.585},

{x1: 31, y1: 0.395, x2: 31, y2: 0.605},

{x1: 32, y1: 0.4115, x2: 32, y2: 0.5885},

{x1: 33, y1: 0.3905, x2: 33, y2: 0.6095},

{x1: 34, y1: 0.409, x2: 34, y2: 0.591},

{x1: 35, y1: 0.4045, x2: 35, y2: 0.5955},

{x1: 36, y1: 0.401, x2: 36, y2: 0.599},

{x1: 37, y1: 0.41, x2: 37, y2: 0.59},

{x1: 38, y1: 0.43, x2: 38, y2: 0.57},

{x1: 39, y1: 0.409, x2: 39, y2: 0.591},

{x1: 40, y1: 0.4235, x2: 40, y2: 0.5765},

{x1: 41, y1: 0.416, x2: 41, y2: 0.584},

{x1: 42, y1: 0.434, x2: 42, y2: 0.566},

{x1: 43, y1: 0.434, x2: 43, y2: 0.566},

{x1: 44, y1: 0.427, x2: 44, y2: 0.573},

{x1: 45, y1: 0.4335, x2: 45, y2: 0.5665},

{x1: 46, y1: 0.4265, x2: 46, y2: 0.5735},

{x1: 47, y1: 0.4335, x2: 47, y2: 0.5665},

{x1: 48, y1: 0.4215, x2: 48, y2: 0.5785},

{x1: 49, y1: 0.4155, x2: 49, y2: 0.5845},

{x1: 50, y1: 0.446, x2: 50, y2: 0.554},

{x1: 51, y1: 0.437, x2: 51, y2: 0.563},

{x1: 52, y1: 0.451, x2: 52, y2: 0.549},

{x1: 53, y1: 0.4515, x2: 53, y2: 0.5485},

{x1: 54, y1: 0.46, x2: 54, y2: 0.54},

{x1: 55, y1: 0.45, x2: 55, y2: 0.55},

{x1: 56, y1: 0.43, x2: 56, y2: 0.57},

{x1: 57, y1: 0.4545, x2: 57, y2: 0.5455},

{x1: 58, y1: 0.458, x2: 58, y2: 0.542},

{x1: 59, y1: 0.466, x2: 59, y2: 0.534},

{x1: 60, y1: 0.4645, x2: 60, y2: 0.5355},

{x1: 61, y1: 0.466, x2: 61, y2: 0.534},

{x1: 62, y1: 0.4655, x2: 62, y2: 0.5345},

{x1: 63, y1: 0.468, x2: 63, y2: 0.532},

{x1: 64, y1: 0.472, x2: 64, y2: 0.528},

{x1: 65, y1: 0.468, x2: 65, y2: 0.532},

{x1: 66, y1: 0.483, x2: 66, y2: 0.517},

{x1: 67, y1: 0.4875, x2: 67, y2: 0.5125},

{x1: 68, y1: 0.499, x2: 68, y2: 0.501},

{x1: 69, y1: 0.4875, x2: 69, y2: 0.5125},

{x1: 70, y1: 0.505, x2: 70, y2: 0.495},

{x1: 71, y1: 0.496, x2: 71, y2: 0.504},

{x1: 72, y1: 0.5115, x2: 72, y2: 0.4885},

{x1: 73, y1: 0.523, x2: 73, y2: 0.477},

{x1: 74, y1: 0.5255, x2: 74, y2: 0.4745},

{x1: 75, y1: 0.517, x2: 75, y2: 0.483},

{x1: 76, y1: 0.533, x2: 76, y2: 0.467},

{x1: 77, y1: 0.548, x2: 77, y2: 0.452},

{x1: 78, y1: 0.551, x2: 78, y2: 0.449},

{x1: 79, y1: 0.532, x2: 79, y2: 0.468},

{x1: 80, y1: 0.555, x2: 80, y2: 0.445},

{x1: 81, y1: 0.5555, x2: 81, y2: 0.4445},

{x1: 82, y1: 0.583, x2: 82, y2: 0.417},

{x1: 83, y1: 0.597, x2: 83, y2: 0.403},

{x1: 84, y1: 0.5935, x2: 84, y2: 0.4065},

{x1: 85, y1: 0.629, x2: 85, y2: 0.371},

{x1: 86, y1: 0.627, x2: 86, y2: 0.373},

{x1: 87, y1: 0.6495, x2: 87, y2: 0.3505},

{x1: 88, y1: 0.663, x2: 88, y2: 0.337},

{x1: 89, y1: 0.7035, x2: 89, y2: 0.2965},

{x1: 90, y1: 0.702, x2: 90, y2: 0.298},

{x1: 91, y1: 0.743, x2: 91, y2: 0.257},

{x1: 92, y1: 0.785, x2: 92, y2: 0.215},

{x1: 93, y1: 0.819, x2: 93, y2: 0.181},

{x1: 94, y1: 0.8575, x2: 94, y2: 0.1425},

{x1: 95, y1: 0.894, x2: 95, y2: 0.106},

{x1: 96, y1: 0.947, x2: 96, y2: 0.053},

{x1: 97, y1: 0.982, x2: 97, y2: 0.018},

{x1: 98, y1: 0.992, x2: 98, y2: 0.008},

{x1: 99, y1: 0.999, x2: 99, y2: 0.001},

{x1: 100, y1: 1, x2: 100, y2: 0}]