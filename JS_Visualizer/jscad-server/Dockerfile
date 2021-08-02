FROM node:6

# get the openjscad source
RUN git clone https://github.com/jscad/OpenJSCAD.org /openjscad

# change to the corresponding directory
WORKDIR /openjscad

# install the software
RUN npm install

# prepare the bootstrap
RUN npm run bootstrap

ADD run.js /openjscad/packages/cli/run.js

CMD node /openjscad/packages/cli/run.js