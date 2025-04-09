// Add syncronous scolling events for the alginment table.
var isSyncingLeftScroll = false;
var isSyncingRightScroll = false;
var leftDiv = document.getElementById('bothNameTables');
var rightDiv = document.getElementById('bothAlignmentTables');

leftDiv.onscroll = function() {
  if (!isSyncingLeftScroll) {
    isSyncingRightScroll = true;
    rightDiv.scrollTop = this.scrollTop;
  }
  isSyncingLeftScroll = false;
}

rightDiv.onscroll = function() {
  if (!isSyncingRightScroll) {
    isSyncingLeftScroll = true;
    leftDiv.scrollTop = this.scrollTop;
  }
  isSyncingRightScroll = false;
}

// place bars at signaturebar that mark position of binary, asym., noisy
function placeBars(coordinates, diversity_list) {
  // if the current index is binary, asymetric or noisy. do nothing if not.
  // if it it something of these three, then place a small bar in the desired
  // color at the correct position at the signature bar.
  // position = x coordinates of point + offset of canvas
  // (coordinates of points are relative to the canvas, not the screen)

  // i,j,z are indexes of binary,asymetric,noisy array
  var i = 0;
  var j = 0;
  var z = 0;
  var y = document.getElementById('signaturebar').offsetTop;

  for(var x = 1; x <= coordinates.length; x++) {
      if(i < diversity_list[0].length && diversity_list[0][i] == x) {
        var bar = document.createElement("div");
        bar.className = "shannonBar";
        bar.style.background = bin_color;
        bar.style.left = coordinates[x-1].x + document.getElementById("graph").offsetLeft + "px";
        bar.style.top = y + "px";
        bar.id = "binary " + i;
        document.getElementById("bar-graph").append(bar);

        var binary_marker = document.createElement("div");
        binary_marker.className = "binaryMarker";
        binary_marker.style.left = (coordinates[x-1].x - 3.0) + document.getElementById("graph").offsetLeft + "px";
        binary_marker.style.top = (y - 12) + "px";
        binary_marker.id = "binary_marker " + i;
        binary_marker.title = x + " (" + nucleotide_list[0][i][1] + ")";
        i++;
        document.getElementById("bar-graph").append(binary_marker);
      }
      if(j < diversity_list[1].length && diversity_list[1][j] == x) {
        var bar = document.createElement("div");
        bar.className = "shannonBar";
        bar.style.background = asym_color;
        bar.style.left = coordinates[x-1].x + document.getElementById("graph").offsetLeft +"px";
        bar.style.top = y + "px";
        bar.id = "asymmetric " + j;
        document.getElementById("bar-graph").append(bar);

        var asymmetric_marker = document.createElement("div");
        asymmetric_marker.className = "asymmetricMarker";
        asymmetric_marker.style.left = (coordinates[x-1].x - 3.0) + document.getElementById("graph").offsetLeft + "px";
        asymmetric_marker.style.top = (y - 12) + "px";
        asymmetric_marker.id = "asymmetric_marker " + j;
        asymmetric_marker.title = x + " (" + nucleotide_list[1][j][1] + ")";
        j++;
        document.getElementById("bar-graph").append(asymmetric_marker);
      }
      if(z < diversity_list[2].length && diversity_list[2][z] == x) {
        var bar = document.createElement("div");
        bar.className = "shannonBar";
        bar.style.background = noisy_color;
        bar.style.left = coordinates[x-1].x + document.getElementById("graph").offsetLeft +"px";
        bar.style.top = y + "px";
        bar.id = "noisy " + z;
        z++;
        document.getElementById("bar-graph").append(bar);
      }
  }
}

// updates bars at signaturebar that mark position of binary, asym., noisy
function updateBars(coordinates, diversity_list) {
  // if the current index is binary, asymetric or noisy. do nothing if not.
  // if it it something of these three, then place a small bar in the desired
  // color at the correct position at the signature bar.
  // position = x coordinates of point + offset of canvas
  // (coordinates of points are relative to the canvas, not the screen)

  // i,j,z are indexes of binary,asymetric,noisy array
  var i = 0;
  var j = 0;
  var z = 0;
  var y = document.getElementById('signaturebar').offsetTop;

  for(var x = 0; x < coordinates.length; x++) {
      if(i < diversity_list[0].length && diversity_list[0][i] == x) {
        var bar = document.getElementById('binary ' + i);
        var binary_marker = document.getElementById("binary_marker " + i);
        bar.style.left = coordinates[x].x + document.getElementById("graph").offsetLeft + "px";
        bar.style.top = y + "px";
        binary_marker.style.left = (coordinates[x].x + document.getElementById("graph").offsetLeft - 2.5) + "px";
        i++;
      }
      if(j < diversity_list[1].length && diversity_list[1][j] == x) {
        var bar = document.getElementById('asymmetric ' + j);
        var asymmetric_marker = document.getElementById("asymmetric_marker " + j);
        bar.style.left = coordinates[x].x + document.getElementById("graph").offsetLeft + "px";
        bar.style.top = y + "px";
        asymmetric_marker.style.left = (coordinates[x].x + document.getElementById("graph").offsetLeft - 2.5) + "px";
        j++;
      }
      if(z < diversity_list[2].length && diversity_list[2][z] == x) {
        var bar = document.getElementById('noisy ' + z);
        bar.style.left = coordinates[x].x + document.getElementById("graph").offsetLeft + "px";
        bar.style.top = y + "px";
        z++;
      }
  }
}

// color index cells according to type (binary = green, asym.=orange, noisy=red)
function colorChars(diversity_list) {
  for(var i = 1; i <= alignment_query_length; i++) {

    var str = document.getElementById(i + "B").innerText;
    var new_html = "";
    var x = 0;
    var y = 0;
    var z = 0;

    for(var j = 1; j <= index_list.length; j++) {
      if(x < diversity_list[0].length && j == diversity_list[0][x]) {
        new_html += "<span title=\"" + diversity_list[0][x] + "\" style=\"background-color: " + bin_color + "\">" + str.charAt(j-1) + "</span>";
        x++;
      }
      else if(y < diversity_list[1].length && j == diversity_list[1][y]) {
        new_html += "<span title=\"" + diversity_list[1][y] + "\" style=\"background-color: " + asym_color + "\">" + str.charAt(j-1) + "</span>";
        y++;
      }
      else if(z < diversity_list[2].length && j == diversity_list[2][z]) {
        new_html += "<span title=\"" + diversity_list[2][z] + "\" style=\"background-color: " + noisy_color + "\">" + str.charAt(j-1) + "</span>"
        z++;
      }
      else {
        new_html += str.charAt(j-1);
      }
    }
    document.getElementById(i + "B").innerHTML = new_html;
  }

  for(var i = 1; i <= alignment_reference_length; i++) {

    var str = document.getElementById(i + "A").innerText;
    var new_html = "";
    var x = 0;
    var y = 0;
    var z = 0;

    for(var j = 1; j <= index_list.length; j++) {
      if(x < diversity_list[0].length && j == diversity_list[0][x]) {
        new_html += "<span title=\"" + diversity_list[0][x] + "\" style=\"background-color: " + bin_color + "\">" + str.charAt(j-1) + "</span>";
        x++;
      }
      else if(y < diversity_list[1].length && j == diversity_list[1][y]) {
        new_html += "<span title=\"" + diversity_list[1][y] + "\" style=\"background-color: " + asym_color + "\">" + str.charAt(j-1) + "</span>";
        y++;
      }
      else if(z < diversity_list[2].length && j == diversity_list[2][z]) {
        new_html += "<span title=\"" + diversity_list[2][z] + "\" style=\"background-color: " + noisy_color + "\">" + str.charAt(j-1) + "</span>"
        z++;
      }
      else {
        new_html += str.charAt(j-1);
      }
    }
    document.getElementById(i + "A").innerHTML = new_html;
  }
}
function placeIndexMarkers() {
  var new_html = "";
  var index = 0;
  for(var i = 1; i <= index_list.length; i++) {
    if(i % 10 == 0) {
      new_html += "&darr;";
      index = index + 10;
      index_str = index.toString();
      for(var j = 0; j < index_str.length; j++) {
        new_html += index_str.charAt(j);
        i++;
      }
    }
    else {
      new_html += "&nbsp"; // &nbsp;
    }
  }
  document.getElementById("0A").innerHTML = new_html;
  document.getElementById("0B").innerHTML = new_html;
}

// calculate width and left of signature bar
function updateSignaturebar(coordinates, index_list) {
  // point 0 coordinates are relative to canvas, so point coordinates + offset of canvas = x axis coordinates (where signaturebar begins)
  document.getElementById('signaturebar').style.left = coordinates[0].x + document.getElementById("graph").offsetLeft + "px";

  // width = last point coordinates - first point coordinates
  document.getElementById('signaturebar').style.width = coordinates[index_list.length - 1].x - coordinates[0].x + 'px';
}

function updateGraph() {
  document.getElementById('graph').style.marginTop = document.getElementById('signaturebar').offsetHeight  + 16 + 'px';
}

function filter(str, diversity_list) {
  var length = 0;
  var marker = false;
  var filter_index = -1;

  if(str == 'binary') {
    length = diversity_list[0].length;
    marker = true;
    filter_index = 0;
  }
  else if(str == 'asymmetric') {
    length = diversity_list[1].length;
    marker = true;
    filter_index = 1;
  }
  else if(str == 'noisy') {
    length = diversity_list[2].length;
    filter_index = 2;
  }

  if(filters[filter_index] == false) {
    for(var i = 0; i < length; i++) {
      document.getElementById(str + " " + i).style.display = "none";
      if(marker == true) {
        document.getElementById(str + "_marker " + i).style.display = "none";
      }
    }
    filters[filter_index] = true;
    document.getElementById("filter" + filter_index).style.textDecoration = "line-through";
  }

  else if (filters[filter_index] == true) {
    for(var i = 0; i < length; i++) {
      document.getElementById(str + " " + i).style.display = "block";
      if(marker == true) {
        document.getElementById(str + "_marker " + i).style.display = "block";
      }
    }
    filters[filter_index] = false;
    document.getElementById("filter" + filter_index).style.textDecoration = "none";
  }
}

function filter_chartjs(index) {

  var bool = chart.data.datasets[index].hidden;
  chart.data.datasets[index].hidden = !chart.data.datasets[index].hidden;

  if(bool == true) {
    document.getElementById("filter" + index + "_chartjs").style.textDecoration = "none";
  }
  else {
    document.getElementById("filter" + index + "_chartjs").style.textDecoration = "line-through";
  }
  chart.update();
}

// *** chart of shannon entropy ***
var ctx = document.getElementById('graph').getContext('2d');

var chart = new Chart(ctx, {
    // The type of chart we want to create
    type: 'line',

    // The data for our dataset
    data: {
      datasets: [{
            label: 'Moving Average',
            data: average_list,
            borderColor: 'rgb(0, 0, 0)',
            hidden: false,

            // dont display points, only the line
            fill: false,
            pointRadius: 0,
            backgroundColor: 'rgba(0, 0, 0, 0.1)'
      },{
            label: 'Entropy',
            data: entropy_list,
            borderColor: 'rgb(30,144,255)',
            hidden: false,

            // dont display points, only the line
            fill: false,
            pointRadius: 0,
            backgroundColor: 'rgba(0, 0, 0, 0.1)'
      }],
      labels: index_list
    },

    // Configuration options go here
    options: {
      tooltips: {
        callbacks: {
          title: function(tooltipItem, data) {
            return tooltipItem[0].index+1;
          },
          label: function(tooltipItem, chart) {
            var val= chart.datasets[tooltipItem.datasetIndex].data[tooltipItem.index];
            return val.toFixed(2);
          }
        }
      },
      animations: false,
      legend: {
        display: false
      },

      // chart length/width = 8 to keep the chat small
      aspectRatio: 8,
      scales: {
        
        x:
          {
              position: 'top',
              id: 'primaryX',
              grid: {
                //TODO: figure out how to seperate grid step size from ticks
                display: false
              },
              type: 'category',
              ticks: {
                  autoSkip: false,
                  min: 0,
                  max: index_list.length,
                  // dont show any label
                  callback: function(value, index, values) {
                    if (value % 50 === 0) {
                       return value;
                     } else {
                       return ' ';
                     }
                  },

                  // dont rotate rabels
                  maxRotation: 0,
                  minRotation: 0
              }
            }
            // // second x-axis to display the grid in stepsize 100
            // {
              // position: 'bottom',
              // id: 'secondaryX',
              // type: 'linear',
              // ticks: {
                // min: 1,
                // max: index_list.length,

                // // dont rotate the chart labels (default is 45deg rotated)
                // maxRotation: 0,
                // minRotation: 0
              // }
            // }
         
        }
    }
  });

  // color index cells where index is binary, asym. or noisy
  colorChars(diversity_list);

  setTimeout(function() {
    // get data of points, e.g. x,y coordinates
    var coordinates = chart.getDatasetMeta(0).data;

    // place bars to mark binary,asym. and noisy positions
    placeBars(coordinates, diversity_list);

    // place index markers every 10 steps
    placeIndexMarkers();

    // calculate signaturebar left and top (bar already exists, only needs update)
    updateSignaturebar(coordinates, index_list);
  }, 400);

  // dynamically resize bar when user zooms
  window.addEventListener("resize", function() {
    // wait until graph resized before the signaturebar resize (900ms is enough)
    setTimeout(function(){
        // resize signaturebar
        var coordinates = chart.getDatasetMeta(0).data;
        updateSignaturebar(coordinates, index_list);
        updateBars(coordinates, diversity_list);
        //updateGraph();
    }, 1800);
  });
