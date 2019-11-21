from fastapi import FastAPI, Body
from typing import Dict
from pydantic import BaseModel, Schema

class PlotInput(BaseModel):
    start_date: str = Schema('2019-08-01', description='Set a start date')
    end_date: str = Schema('2019-10-25', description='Set an end date')
    start_hour: str = Schema('6:00', description='Set a starting hour')
    end_hour: str = Schema('11:00', description='Set an ending hour')
    keys_to_compare: list = Schema(..., description='Give a list of keys you want to plot')

example_plot_input = Body(
    ...,
    example={
        'start_date': '2019-08-01',
        'end_date': '2019-10-25',
        'start_hour': '0:00',
        'end_hour': '0:00',
        'keys_to_compare': ['veh', 'TEMP', 'no2']
    }
)