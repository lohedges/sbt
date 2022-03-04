/*
  Copyright (c) 2021-2022 Lester Hedges <lester.hedges@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <vector>

#include <poplar/DeviceManager.hpp>
#include <poplar/IPUModel.hpp>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

/// Setup and connect to a Poplar IPU device.
poplar::Device setIpuDevice();

int main()
{
    // Read events from file.
    EventReader event_reader(std::filesystem::path("data"), true);
    event_reader.readEvents(true);

    // Setup the IPU device handle.
    auto device = setIpuDevice();

    // Setup the graph program.
    std::cout << "Creating graph program...\n";
    SearchByTripletIPU search_by_triplet(std::move(device),
                                         event_reader.getEvents());
                                         

    search_by_triplet.execute();

    std::cout << "Finished!\n";

    return 0;
}

poplar::Device setIpuDevice()
{
    auto dm = poplar::DeviceManager::createDeviceManager();
    auto hwDevices = dm.getDevices(poplar::TargetType::IPU, 4);
    if (hwDevices.size() > 0)
    {
        for (auto &d : hwDevices)
        {
            if (d.attach())
            {
                return std::move(d);
            }
        }
    }

    std::cerr << "Unable to connect to a device with 4 IPUs.\n";
    exit(-1);
}
