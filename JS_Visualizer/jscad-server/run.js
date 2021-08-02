// /openjscad/packages/cli

const fs = require('fs');
const http = require('http');
const generateOutputData = require('./generateOutputData');

let inputFormat = "jscad";
let outputFormat = "stl";
let addMetaData = false;
let version = "1.10.0";
let params = {};
//let src = fs.readFileSync(inputFile, inputFile.match(/\.stl$/i) ? 'binary' : 'UTF8');

//create a server object:
http.createServer(function (request, res) {
    let body = '';
    request.on('data', (chunk) => {
        body += chunk;
    }).on('end', () => {
        generateOutputData(body, params, {undefined, outputFormat, undefined, inputFormat, version, addMetaData})
        .then(function (outputData) {
            res.setHeader("Access-Control-Allow-Origin", "*");
            res.setHeader("Access-Control-Allow-Headers", "X-Requested-With");
            res.write(outputData.asBuffer()); //write a response to the client
            res.end(); //end the response
        }).catch(error => console.error(error));
    });
}).listen(8096); //the server object listens on port 8080