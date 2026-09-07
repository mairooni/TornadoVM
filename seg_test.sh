#!/bin/bash

while true
do
  # Replace "your_command" with the actual command you want to execute
  tornado --jvm "-Xmx6g -Dtornado.recover.bailout=False -Dtornado.unittests.verbose=True -Dtornado.device.memory=2GB"  -m  tornado.unittests/uk.ac.manchester.tornado.unittests.tools.TornadoTestRunner  --params "uk.ac.manchester.tornado.unittests.multithreaded.TestMultiThreadedExecutionPlans"


  # Optionally, add a sleep timer to avoid overwhelming the system
  # sleep 1  # Sleep for 1 second before the next iteration
done
