/* This code is designed to be run on an Arduino board to receive and process the 4x4 sensor array data.. */

#define Receive_buffer_length 16*2
#define Delimiter_size 2
char dCar[Receive_buffer_length + Delimiter_size];

void setup() 
{
  Serial.begin(500000);
  Serial1.begin(500000);
  while (!Serial) {}
  while (!Serial1) {}
  Serial.setTimeout(100);
  Serial1.setTimeout(100);
}

uint8_t _buffer[Receive_buffer_length]={0};

void loop() 
{
    if(Serial1.available()>0)
    {
      Serial1.readBytes(_buffer, Receive_buffer_length+Delimiter_size);
      Serial1.flush();
      for(int i=0; i< Receive_buffer_length; i=i+2){
        Serial.print((_buffer[i+1] << 8) | _buffer[i]);
        Serial.print(' ');
      }
      Serial.println("");
    }
}
